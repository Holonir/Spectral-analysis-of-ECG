Исследование разложения сигнала по методу быстрого преобразования Фурье (БПФ), изучение спектров сигналов ЭКГ в норме и при различных патологиях.
Состав проекта: 
•	Программа, позволющая рассчитать спектральные оценки для фрагментов ЭКГ, соответствующих трем различным видам сердечного ритма.
При исследовании медико-биологических сигналов часто используются методы спектрального анализа, позволяющие получить численные оценки частотного состава сигнала. Наиболее распространен спектральный анализ, основанный на дискретном преобразовании Фурье.
Рассчитанные с использованием БПФ спектральные оценки имеют два существенных недостатка: 
1) в них может присутствовать сильная составляющая на нулевой частоте;
2) из-за необходимости предположения о периодичности сигнала могут возникнуть ложные высокочастотные составляющие, вызванные скачками на краях анализируемого фрагмента.
Для устранения этих нежелательных эффектов из сигнала перед БПФ вычитают среднее значение и далее умножают сигнал на так называемую “оконную” функцию, подавляющую скачки на краях.
