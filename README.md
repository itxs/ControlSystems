# ControlSystems

Набор кода SciLab и схем XCos для решения типовых задач ТАУ.
В дальнейшем может появиться код и схемы Matlab/Simulink, SimInTech, и тд.

### Краткая справка:
#### ProcApprox.sce
Библиотека с исходным кодом SciLab для идентификации систем 1 и 2 порядка с задержкой (FOPDT, SOPDT).
Используется стандартная функция datafit, критерий оценки качества - квадратичный
Для использования подготовьте массив калибровочных измерений с датчика,
подав на управляемый объект ступенчатое воздействие с максимально возможной амплитудой.
Например, если вы управляете мотором, регуляируя ток, подайте на него максимально возможный ток,
и замерьте данные с датчика оборотов или энкодера с требуемым интервалом между измерениями.
Получится массив показаний с  датчика с зависимостью от времени. Загрузите этот массив в SciLab
и создайте там линейный массив точек времени (например с помощью linspace).
Для идентификации объекта управления (мотор) вызовите функцию
```
[coeffs, error] = ProcessApproximate(x, y, "sopdt");
```
Где x - линейный массив меток времени, y - массив данных с датчика.

На выходе получим две переменные, в одной из которых будет вектор коэффициентов
передаточной функции записанной в стандартном виде, а в другой - полученная ошибка сходимости.
Можно выбрать одну из двух стандартных передаточных функций:
```
exp(-b*s)/(a*s^2+c*s+d)
exp(-b*s)/(a*s+c)
```
Есть возможность задать свои начальные значения вектора коэффициентов, если
при стандартных значениях (все к-ты равны 5) сходимости нет. Кроме того, есть возможность задания вектора весовых коэффициентов, если некоторая область
значений с датчика более достоверна чем другая.

Пример использования всех параметров:
```
[coeffs, error] = ProcessApproximate(x, y, "sopdt", [начальные значения коэффициентов], [веса показаний датчика]);
```
