# -*- coding: utf-8 -*-
'''
Created on 13 февр. 2018 г.

Программа для проекта "Спектрометр". Создаёт базу данных *.asc объектов, взятых с сайта http://speclab.cr.usgs.gov.
Для работы необходимо, чтобы в папке запуска находилась папка с объектами имеющая имя "objects" и папка с данными фильтров,
имеющая имя "filters". Используются 10 фильтров комплекта Thorlabs FKB-VIS-40 от 400 до 850 нм.
Файл базы данных в каталоге запуска создаётся автоматически после завершения работы программы и имеет имя "data.csv".
База данных имеет вид следующей таблицы, где A, B - значения интеграла произведения функций фильтра и объекта.
 _____________________________________________________
|             |                |                |     \
|             | Filter_1(400)  | Filter_2(450)  | Filte\
|_____________|________________|________________|_______\
|             |                |                |       /
| Object_name |       A        |       B        |      |
|_____________|________________|________________|_____/

@author: Greg
'''
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import pandas as pd
import os, collections

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Принимает путь и имя файла с данными. Возвращает словарь. Ключ - имя объекта (15я строка в файле), данные - pd.DataFrame с тремя столбцами ('wavelength', 'reflectance', 'deviation').
# При создании pd.DataFrame отсеиваются значения меньше нуля, те, которые распологаются за пределами интервала [361 - 890 нм] и значения имеющие отриццательный коэфф. отражения (см. файл с данными)
# Этот интервал - макс. и мин. значения длин фолн фильтров. Существует два вида файлов с данными, которые отличаются представлением шкалы длин волн.
# Функция это учитывает и автоматически определяет множитель (100 или 1000) для перевода шкалы в нм.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
def obj_dataset(path, filename):
    with open(path + '\\' + filename) as file:
        data = file.readlines()
        
    dataset = pd.DataFrame(columns=('wavelength', 'reflectance', 'deviation'))
    WavelengthFactorFlag = False # Флаг показывает определён множитель шкалы длин волн или нет. Множитель определяется по первому неотрицательному значению длины волны.
    
    for i in range(16, len(data)):
        datastring = data[i].strip().split('      ') # строка с данными из файла объекта
        
        if float(datastring[0]) < 0: # отсеиваем отрицательные длины волн
            continue
            
        if not WavelengthFactorFlag: # определяемем множитель шкалы длин волн
            WavelengthFactor = 1000 if(int(datastring[0][0]) == 0) else 100
            WavelengthFactorFlag = True
                    
        if ((float(datastring[0]) * WavelengthFactor) >= 361 and ((float(datastring[0]) * WavelengthFactor) <= 890) and float(datastring[1]) > 0): # отсеиваем строки вне интервала и со значением < 0
            
            dataset = dataset.append(pd.DataFrame([[int(float(datastring[0]) * WavelengthFactor),
                                                    float(datastring[1]),
                                                    float(datastring[2])]],
                                                    columns=('wavelength', 'reflectance', 'deviation')),
                                                    ignore_index=True)
        else:
            continue
    return {data[14].strip() : dataset}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Принимает путь до папки с файлами фильтров. Возвращает словарь ключи которого - имена фильтров, а данные - pd.DataFrame с двумя столбцами ('wavelength', 'transmission').
# Используется упорядоченная версия словаря из collections
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
def filter_database(path):
    fltDataBase = collections.OrderedDict()
    files = sorted([file for file in os.listdir(path) if file.lower().endswith('.xlsx')], key=len)
    for file in files:
        get_fil = pd.read_excel(path + '\\' + file)
        filx = get_fil['Wavelength (nm)']
        fily = get_fil['% Transmission']
        dataset = pd.concat([filx, fily], axis=1, ignore_index=True)
        dataset.columns = ['wavelength', 'transmission']
        fltDataBase.update({file.split('.')[0] : dataset})
    return fltDataBase
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Принимает: flt - значения упорядоченного словаря, являющиеся pd.DataFrame с данными итерационно выбранного фильтра.
# obj - значения словаря, являющиеся pd.DataFrame с данными объекта (в словаре объекта всегда один элемент).
# Возвращает значение интеграла результирующей функции.
# Данные фильтра интреполируются. В полученную функцию подставляются обрезанные до интервала [мин. длины волны фильтра - макс. длины волны фильтра]
# данные образца для получения данных фильтра в такой же сетке (по длинам волн). Далее значения пропускания фильтра и отражения образца перемножаются
# и полученная функция интегрируется.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
def get_integral_value(flt, obj):

    fltFunction = interpolate.interp1d(flt.wavelength, flt.transmission) # Интерполируем данные фильтра

    fltDataMin = min(flt.wavelength) # Находим макс. и мин. длины волн фильтра
    fltDataMax = max(flt.wavelength)

    cropObj_data = obj[(obj.wavelength >= fltDataMin) & (obj.wavelength <= fltDataMax)].reset_index(drop=True) # Обрезаем данные объекта по мин. и макс. полученным выше
    
    newxFltData = cropObj_data.wavelength.reset_index(drop=True) # Сбрасываем индексы (чтобы шли с 1)
    newxFltData = pd.to_numeric(newxFltData) # Преобразуем последовательности obj  в float64 <--- Возможный костыль! Мб можно сразу создать норм. DataFrame 
    
    newyFltData = fltFunction(newxFltData) # Получаем значения прохождения для фильтра в такойже сетке как и у объекта
    
    newFlt_data = pd.DataFrame({'wavelength' : newxFltData, 'transmission' : newyFltData}) # Создаём pd.DataFrame с новыми данными для фильтра
  
    prodFltObj = pd.DataFrame({'wavelength' : newFlt_data.wavelength, 'reflectance' : newFlt_data.transmission * cropObj_data['reflectance']}) # Перемножаем значения пропускания фильтра и отражения образца
    
    result = integrate.trapz(prodFltObj.reflectance, x = prodFltObj.wavelength) # Интегрируем полученные данные
    
    return result
#-------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    PATH = os.getcwd()
    #PATH = 'C:\\dzz\\splib06.library'
    FLT_PATH = PATH + '\\' + 'filters'
    OBJ_PATH = PATH + '\\' + 'objects'

    filters = filter_database(FLT_PATH)  # упорядоченный словарь с всеми фильтрами
    ObjDataBase = pd.DataFrame()  # pd.DataFrame для записи результата работы программы
    objFiles = sorted([file for file in os.listdir(OBJ_PATH) if file.lower().endswith('.asc')], key=len) # Список объектов в рабочей папке
    
    for file in objFiles:
        
        os.system('cls')
        print(u'Выполнено %d' % (objFiles.index(file) * 100 / len(objFiles)) + u'%', end = '', flush = True)
        
        objDict = obj_dataset(OBJ_PATH, file) # Получаем словарь объекта
        objName = list(objDict.keys())[0]     # Имя объекта
        objData = list(objDict.values())[0]   # Данные объекта

        ObjIntegralValue = []
        for flt in filters:
            ObjIntegralValue.append("{0:.3f}".format(get_integral_value(filters.get(flt), objData)))  # Получаем значение интеграла, обрезаем до .3f и добавляем в ObjIntegralValue.
        ObjDataBase = ObjDataBase.append(pd.DataFrame([ObjIntegralValue], columns=['Filter_1(400)', 'Filter_2(450)', 'Filter_3(500)', 'Filter_4(550)', 'Filter_5(600)',
                                                                               'Filter_6(650)', 'Filter_7(700)', 'Filter_8(750)', 'Filter_9(800)', 'Filter_10(850)'],
                                                                               index=[objName]))
ObjDataBase.to_csv('data.csv') # сохраняем в рабочей папке
#-------------------------------------------------------------------------------------------------------

    
