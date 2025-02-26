# Программа расчёта координат для 20 молекул 
## Введение 
Изучение взаимодействия молекул играет важную роль в объяснении свойств различных веществ и их состояний. Наблюдение за поведением молекул позволяет определять макроскопические характеристики субстанций и строить модели состояния веществ в определённых условиях.   
В данной работе представлена модель движения 20 взаимодействующих  молекул. 
В работе используется библиотека SFML Graphics.

## 1. Цели и задачи 

Целью данной работы является создание программной реализации модели взаимодействия молекул, вычисляющей  изменения положений молекул и силы влияния их друг на друга. Для её достижения были поставлены следующие задачи: 

1) Ознакомиться с основными силами взаимодействия молекул и способами вычисления этих сил. 
2) Научиться вычислять величину результирующего влияния  молекул друг на друга. 
3) Изучить способы вычисления изменения положения молекул вследствие их взаимодействия. 
4) Научиться определять расположение молекул в пространстве в определённый момент времени . 
5) Построить модель взаимодействия молекул для решения выбранной задачи. 
6) Провести несколько экспериментов на различных наборах входных данных. 
7) Сделать программную реализацию. 


## 2. Обзор


### 2.1 Природа взаимодействия молекул 

Межмолекулярное взаимодействие имеет электростатическую природу. Предположение о его существовании было впервые использовано Йоханнесом Ван дер Ваальсом в 1873 году для объяснения свойств реальных и жидкостей. В широком смысле под ним можно понимать взаимодействия между любыми частицами (молекулами, атомами), при которых не происходит образования химических связей. Ранее в 1869 году   Ван дер Ваальс открыл сами силы межмолекулярного взаимодействия, позже названные Вандерваальсовыми силами. 
Модель построенная в данной работе основывается на Потенциале  Леннард-Джонса, описывающем парное взаимодействие неполярных молекул, а также зависимость энергии взаимодействия двух частиц от расстояния между ними. Эта модель достаточно реалистично передаёт свойства реального взаимодействия сферических неполярных молекул и поэтому широко используется в расчётах и при компьютерном моделировании. 
Молекулы попарно взаимодействуют друг с другом, для нахождения значения и направления ускорения каждой молекулы необходимо вычислять результирующую силу взаимодействия между ней и остальными рассматриваемыми молекулами. 


### 2.2 ’’Вылет’’ молекул за пределы радиуса взаимодействия 

В модели потенциала Леннарда-Джонса молекулы взаимодействуют  таким образом, что при увеличении расстояния между двумя молекулами  их ускорения меняют направление в сторону друг друга, тем самым оба тела начинают сближаться. Аналогично с уменьшением расстояния: чем ближе молекулы, тем больше сила, с которой они отталкиваются друг от друга.  

Однако при сильном сближении молекул сила Ван-Дер-Ваальса придаёт им гораздо более сильное ускорение, чем при большом  отдалении друг от друга, что приводит к вылету некоторых молекул из за пределы эффективного радиуса взаимодействия. 


### 2.3 Язык программирования 

Язык С++ является расширением языка С, позволяющий использовать такие инструменты как объектно-ориентированное программирование и шаблоны. Он  удобен и имеет широкий спектр применения в создании программного обеспечения. Код на С++ быстр и эффективен.  

 

## 3. Постановка задачи 

Для выполнения работы была поставлена задача вычислить координаты столкновения с поверхностью сосуда двадцати молекул газа движущихся с определённой скоростью. 
Как пример газа мною был выбран кислород. Начальное положение молекул заключено в сферу, радиус которой равен одному нанометру. Расстояние от молекул до стенки сосуда равен 10 радиусам взаимодействия. Начальную скорость движения молекул по оси Y задаёт пользователь программы, начальная скорость по остальным осям равна 0. 


## 4. Ход работы 


### 4.1 Начальные условия 

Создаётся класс Molecule, объекты которого обладают такими атрибутами как speedx,  speedy,  speedz, отображающие скорость молекул по трём координатным осям, а также xcoord, ycoord и zcoord, обозначающие положение молекул в пространстве. 
Для построения модели движения молекул необходимо задать начальные значения скорости и координат для каждой молекулы. Пользователь вводит её среднее начальное значение, варьирующееся от 0 до 331 метрам в секунду (последнее значение примерно соответствует  скорости звука), после чего с помощью функции нормального распределения каждая молекула получает свою собственную скорость. 
Начальное положение молекул определяется с помощью функции равномерного распределения в сфере радиуса n на расстоянии 10n до поверхности сосуда, которая соответствует нулевым координатам оси Y. 
Длина n  одной клетки координатной плоскости равна 1 нанометру. 
Начальное положение центра сферы по осям X и Y равно 0. 
Также определяются коэффициенты ε, σ и масса молекулы вещества. 
В случае кислорода ε=1.629166 * 10^(-21) Дж 
σ=3.52 *10^(-10) м 
mas=5.33 * 10^(-26) кг. 

 
### 4.2 Результирующая сила взаимодействия 

Для вычисления результирующей силы влияния на одну молекулу последовательно высчитывается расстояние между выбранной молекулой и одной из оставшихся, после чего передаётся в функцию lennard_jones_force, где рассчитывается сила взаимодействия двух молекул через взятие производной потенциала Леннард-Джонса. Полученное значение передаётся в вектор сил взаимодействия между выбранной молекулой и остальными. После заполнения вектора сил идёт следующий шаг. 

Результирующая сила вычисляется с помощью функции resultforce. 

Силы влияния двух молекул друг на друга представляется как вектор в пространстве. Для того чтобы было возможно сложить векторного массива сил в выбранном координатном пространстве берётся коэффициент k, равный отношению силы взаимодействия двух молекул на расстояние между ними. Далее берётся сумма сил влияния на выбранную молекулы со стороны оставшихся следующим образом: 
1) Берётся разность между координатами двух молекул по одной из координатных осей. 
2) Полученное значения умножается на коэффициент k и прибавляется к сумме сил по выбранной координатной оси. 
3) Эти шаги проделываются для всех координатных осей и каждого элемента векторного массива сил, в итоге получаем значения результирующей силы действия всех молекул на одну выбранную по 3 координатам. 
4) Полученные значения суммируются с начальной силой молекулы, обоснованной её начальной скоростью, и передаются в вектор класса Forcevector, объекты которого обладают тремя атрибутами отвечающими за проекцию сил на координатные оси. Таким образом мы получили результирующую силу для одной молекулы. 
Проделав шаги 1-4 для каждой из молекул получаем заполненный вектор результирующих сил всех молекул. 



### 4.3 Изменение положения молекул 

Вычисление изменений координат молекул: 
1) Рассчитываем проекции ускорения молекулы на координатные оси  по формуле  a=F/m, где F-один из атрибутов объекта вектора результирующих сил. 
2) Меняем значения проекций скорости и координат молекулы по формулам. 
3) Применяем шаги 1-2 для каждой молекулы. 
Полученные значения соответствуют изменению положения молекул за 1 пикосекунду. 

 

### 4.4 Нахождение точек столкновения молекул с поверхностью 

Нужно найти координаты по осям X и Z в момент когда координата по Y будет равна 0. Для этого проверяем значения координаты Y каждой молекулы после изменения положения за пикосекунду. Если значение Y стало отрицательным, последнее изменение атрибутов молекулы откатывается и вычисляется время t1, через которое координата Y станет равна 0, подставив которое получаются координаты столкновения молекулы с поверхностью. 
Столкновение с поверхностью происходит  абсолютно упруго, то есть импульс и внутренняя энергия молекул остаются неизменными, дальнейшее движение молекулы продолжается в противоположном направлении по координате Y, направление по осям X и Z не меняется. 
Программа работает в цикле в течение 1 наносекунды, поскольку за это время каждая молекула либо столкнулась с поверхностью, либо вылетела за пределы взаимодействия. Каждый шаг цикла отображает одну пикосекунду движения молекул. Значения координат столкновения записываются в вектор Final. 
Для того чтобы записывалось только первая точка столкновения молекулы проводится проверка одноразового заполнения вектора Final через вспомогательный векторный массив Check. 
Программа выводит на экран координаты соприкосновения молекул со стенкой сосуда, а также график изменения положения молекул на примере графика движения трёх молекул в проекции на координатную плоскость X0Y. 


### Заключение 

В данной работе были получены следующие результаты: 
1)Изучены основные силами взаимодействия молекул, способы вычисления сил и механика изменения положения молекул в пространстве. 
2)Построена система вычисления положения взаимодействующих молекул. 
3)Проведены эксперименты на нескольких наборах входных данных. 
4) Сделана программная реализация на языке С++. 
