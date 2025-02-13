#include <iostream>
#include <stdio.h>
#include <random>
#include <iomanip>
#include <time.h>
#include <string>
#include <map>
#include <SFML/Graphics.hpp>
#include <cmath>

using namespace std;
const int W = 1080;
const int H = 720;
//класс молекула
class Molecule
{
public:
	long double  speedx;
	long   double  speedy;
	long    double  speedz;
	long double  xcoord;
	long double  ycoord;
	long double  zcoord;
};
//класс вектор силы
class Forcevector
{
public:
	long  double  coordx;
	long  double  coordy;
	long  double  coordz;
};

//начальное распределение скоростей,используется нормальное распределение
void distribution(const double m, const double s, const int samples, vector <Molecule>& pack)
{
	mt19937 gen(1701);
	normal_distribution<> distr(m, s);


	// generate the distribution as a histogram
	map<long double, int> histogram;
	for (int i = 0; i < samples; ++i) {
		++histogram[distr(gen)];
	}
	int counter = 0;
	for (const auto& elem : histogram) {
		pack[counter].speedy = (long double)elem.first;
		pack[counter].speedy = (-1) * pack[counter].speedy;
		counter++;
	}


}
//расстояние между молекулами в координатах
long double  Distance(long double  x, long double  y, long double  z)
{
	long double  r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	return r;

}
//результирующий вектор силы для одной молекулы
void resultforce(vector <long double >& Forcepack, int m, vector <Molecule>& Workpack, vector <Forcevector>& Vectors, long double distconst)
{
	long double  fx = 0;
	long double  fy = 0;
	long double  fz = 0;
	long double  r = 0;
	long double  k = 0;

	for (int i = 0; i < 20; i++)
	{
		r = Distance(Workpack[i].xcoord - Workpack[m].xcoord, Workpack[i].ycoord - Workpack[m].ycoord, Workpack[i].zcoord - Workpack[m].zcoord);
		//k-коэффициент, ставящий в соответсвие величине силы взаимодействия двух молекул расстояние между ними,
		//тем самым мы получаем вектор силы с нужным направлением, а его длина равна самой силе взаимодействия,
		//благодаря чему мы можем брать векторную сумму
		if (r == 0) k = 0;
		else k = Forcepack[i] / r;
		fx = fx + (Workpack[i].xcoord - Workpack[m].xcoord) * k;
		fy = fy + (Workpack[i].ycoord - Workpack[m].ycoord) * k;
		fz = fz + (Workpack[i].zcoord - Workpack[m].zcoord) * k;
	}
	//прибавляем к начальному вектору силы молекулы результирующий вектор взаимодействия с остальными молекулами
	Vectors[m].coordx = Vectors[m].coordx + fx;
	Vectors[m].coordy = Vectors[m].coordy + fy;
	Vectors[m].coordz = Vectors[m].coordz + fz;

}



//сила взаимодействия двух молекул
long double  lennard_jones_force(long double  r2, long double sigma, long double  eps)
{
	long double  inv_r2;
	if (r2 == 0) inv_r2 = 0;
	else inv_r2 = 1 / r2;
	long  double  sigma_over_r6 = pow(sigma * inv_r2, 6);
	long  double  sigma_over_r12 = pow(sigma * inv_r2, 12);
	long  double   force = 24 * eps * inv_r2 * (2 * sigma_over_r12 - sigma_over_r6);
	long double  potential = 4 * eps * (sigma_over_r12 - sigma_over_r6);
	return force;
}


long double  GetRandomNumberFloat(double  min, double max, int precision)
{
	static std::mt19937 gen(time(NULL));
	std::uniform_real_distribution<> ufd(min, max);
	return ufd(gen);
}

//функция графического отображения движения 3 молекул

void Graphics(vector <long double >& Molecule1, vector <long double >& Molecule2, vector <long double >& Molecule3,int count )
{
	sf::RenderWindow window(sf::VideoMode(W, H), "Function graph!");

	int x0 = W / 2;
	int y0 = H / 2;

	sf::CircleShape point(2.f);
	point.setFillColor(sf::Color::Blue);

	int sc = 40;

	sf::RectangleShape line[40];
	for (int i = 0; i < 40; i++) {
		line[i].setSize(sf::Vector2f(1, 20));
		line[i].setFillColor(sf::Color::Black);

		if (i < 20) {
			if (i < 10)
				line[i].setPosition(x0 - (i + 1) * sc, y0 - 10);
			else
				line[i].setPosition(x0 + (i - 9) * sc, y0 - 10);
		}
		else {
			line[i].setRotation(90);
			if (i < 30)
				line[i].setPosition(x0 + 10, y0 + (i - 30) * sc);
			else
				line[i].setPosition(x0 + 10, y0 + (i - 29) * sc);
		}
	}

	sf::RectangleShape OsX(sf::Vector2f(W, 1));
	OsX.setFillColor(sf::Color::Black);
	OsX.setPosition(0, y0);

	sf::RectangleShape OsY(sf::Vector2f(1, H));
	OsY.setFillColor(sf::Color::Black);
	OsY.setPosition(x0, 0);

	sf::RectangleShape strel[4];
	for (int i = 0; i < 4; i++) {
		strel[i].setSize(sf::Vector2f(1, 25));
		strel[i].setFillColor(sf::Color::Black);

		if (i < 2)
			strel[i].setPosition(x0, 0);
		else
			strel[i].setPosition(W, y0);
	}
	strel[0].setRotation(25);
	strel[1].setRotation(-25);
	strel[2].setRotation(60);
	strel[3].setRotation(-250);


	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		
		window.clear(sf::Color::White);
		window.draw(OsX);
		window.draw(OsY);
		for (int i = 0; i < 4; i++)
			window.draw(strel[i]);
		for (int i = 0; i < 40; i++)
			if (i != 19 && i != 20)
				window.draw(line[i]);
		
		for (int i = 0; i < count; i = i + 2) 
		{
			
			float x1 = x0+ Molecule1[i]*10 ;
			float y1 =  y0-Molecule1[i+1]*10 ;
			point.setPosition(x1, y1);
			point.setFillColor(sf::Color::Black);
			window.draw(point);	
		}
		
		for (int i = 0; i < count; i = i + 2)
		{
			float x2 = x0 + Molecule2[i] *10;
			float y2 = y0 - Molecule2[i + 1] *10;
		point.setPosition(x2, y2);
		point.setFillColor(sf::Color::Blue);
		window.draw(point);
		}
		for (int i = 0; i < count; i = i + 2)
		{

			float x3 = x0 + Molecule3[i] *10;
			float y3 = y0 - Molecule3[i + 1] *10;
			point.setPosition(x3, y3);
			point.setFillColor(sf::Color::Red);
			window.draw(point);
		}
		window.display();
	}

	
}
int main()
{
	long double radius = 2.2 * pow(10, -10);// эффективное расстояние взаимодействия молекул
	long double r = 0;
	long double t = 1 * pow(10, -12);
	long double distconst = 1 * pow(10, -9);
	long double  t1 = 0;
	long double  t2 = 0;
	long double  mas = 5.33 * pow(10, -23);//масса молекулы кислорода
	long double  sigma = 0.352 * pow(10, -9);
	long double  eps = 1.629166 * pow(10, -21);//глубина потенциальной ямы
	float s;
	cout << "Consider the movement of 20 oxygen molecules";
	cout << endl;
	cout << "Input basic speed,  meter/second (0-331)";
	cout << endl;
	cin >> s;
	//векторы координат трёх молекул
	vector <long double> Molecule1;
	vector <long double> Molecule2;
	vector <long double> Molecule3;
	vector <Molecule> Workpack(20);// вектор обьектов класса Молекула
	vector <long double > Forcepack(20);// вектор сил действующих на молекулу
	vector <Forcevector> Vectors(20);// вектор результирующих сил молекул
	vector <Molecule> Final(20);// вектор координат столкновения с поверхностью
	vector <int> Check(20);// проверка единоразового заполнения вектора Final
	//распределяем скорости по начальным данным
	distribution(s, s / 4, 20, Workpack);

	//задаём начальное положение молекул в координатах, расстояние в 10 радиусов взаимодействия отображается в координате у
	vector <long double> Initialx(20);
	vector <long double> Initialy(20);
	vector <long double> Initialz(20);
	for (int i = 0; i < 20; i++)
	{
		Workpack[i].xcoord = GetRandomNumberFloat(-1, 1, 10) * distconst;
		Initialx[i] = Workpack[i].xcoord;
		Workpack[i].ycoord = GetRandomNumberFloat(-1, 1, 10) * distconst + 10*radius;
		Initialy[i] = Workpack[i].ycoord;
		Workpack[i].zcoord = GetRandomNumberFloat(-1, 1, 10) * distconst;
		Initialz[i] = Workpack[i].zcoord;
	}
	//отображение начальных координат точек
	Molecule1.push_back(Initialx[0]/ distconst);
	Molecule1.push_back(Initialy[0] / distconst);
	Molecule2.push_back(Initialx[2] / distconst);
	Molecule2.push_back(Initialy[2] / distconst);
	Molecule3.push_back(Initialx[5] / distconst);
	Molecule3.push_back(Initialy[5] / distconst);
	Graphics(Molecule1, Molecule2, Molecule3,2);
	cout << endl;
	cout << endl;
	cout << endl;
	for (int f = 0; f < 1000; f++)
	{
		for (int q = 0; q < 20; q++)
		{
			for (int j = 0; j < 20; j++)
			{
				//определяем расстояние между молекулами и их силу взаимодействия через потенциал Леннарада-Джонса
				r = Distance(Workpack[j].xcoord - Workpack[q].xcoord, Workpack[j].ycoord - Workpack[q].ycoord, Workpack[j].zcoord - Workpack[q].zcoord);
				Forcepack[j] = lennard_jones_force(r, sigma, eps);
			}
			//вычисляем результирующий вектор силы для молекулы
			resultforce(Forcepack, q, Workpack, Vectors, distconst);
		}
		//изменяем скорость и положение молекулы
		for (int i = 0; i < 20; i++)
		{
			Workpack[i].speedx = Workpack[i].speedx + Vectors[i].coordx * t / mas;
			Workpack[i].speedy = Workpack[i].speedy + Vectors[i].coordy * t / mas;
			Workpack[i].speedz = Workpack[i].speedz + Vectors[i].coordz * t / mas;
			Workpack[i].xcoord = Workpack[i].xcoord + Workpack[i].speedx * t + Vectors[i].coordx * t * t / (2 * mas);
			Workpack[i].ycoord = Workpack[i].ycoord + Workpack[i].speedy * t + Vectors[i].coordy * t * t / (2 * mas);
			Workpack[i].zcoord = Workpack[i].zcoord + Workpack[i].speedz * t + Vectors[i].coordz * t * t / (2 * mas);
			

			//проверяем достигла ли какая либо молекула стенки сосуда(стенка сосуда обозначается как y=0 на координатах);
			if (Workpack[i].ycoord < 0) 
			{
				

				Workpack[i].xcoord = Workpack[i].xcoord - (Workpack[i].speedx * t + Vectors[i].coordx * t * t / (2 * mas));
				Workpack[i].ycoord = Workpack[i].ycoord - (Workpack[i].speedy * t + Vectors[i].coordy * t * t / (2 * mas));
				Workpack[i].zcoord = Workpack[i].zcoord - (Workpack[i].speedz * t + Vectors[i].coordz * t * t / (2 * mas));

				Workpack[i].speedx = Workpack[i].speedx - Vectors[i].coordx * t / mas;
				Workpack[i].speedy = Workpack[i].speedy - Vectors[i].coordy * t / mas;
				Workpack[i].speedz = Workpack[i].speedz - Vectors[i].coordz * t / mas;

				//если координата у молекулы стала отрицательной то мы возврщаемся к предыдущему положению молекулы и вычисляем время когда у стало равно 0
				if ((-2 * Workpack[i].speedy - sqrt(pow(2 * Workpack[i].speedy, 2) - 8 * Vectors[i].coordy / mas * Workpack[i].ycoord)) / (2 * Vectors[i].coordy / mas) > 0)
					t1 = (-2 * Workpack[i].speedy - sqrt(pow(2 * Workpack[i].speedy, 2) - 8 * Vectors[i].coordy / mas * Workpack[i].ycoord)) / (2 * Vectors[i].coordy / mas);

				else t1 = (-2 * Workpack[i].speedy + sqrt(2 * pow(Workpack[i].speedy, 2) - 8 * Vectors[i].coordy / mas * Workpack[i].ycoord)) / (2 * Vectors[i].coordy / mas);
				t2 = t - t1;
				t1 = t;
				//записываем координаты столкновения молекулы с поверхностью
				if (Check[i]  == 0) {
				Final[i].xcoord = Workpack[i].xcoord + (Workpack[i].speedx * t1 + Vectors[i].coordx * t1 * t1 / (2 * mas));
				Final[i].ycoord = 0;
				Final[i].zcoord = Workpack[i].zcoord + (Workpack[i].speedz * t1 + Vectors[i].coordz * t1 * t1 / (2 * mas));
				}
				//записываем изменение координат 3 выбранных для примера молекул
				if (i == 0) {
					Molecule1.push_back(Final[i].xcoord / distconst);
					Molecule1.push_back(Final[i].ycoord / distconst);
				}
				if (i == 2) {
					Molecule2.push_back(Final[i].xcoord / distconst);
					Molecule2.push_back(Final[i].ycoord / distconst);
				}
				if (i == 5) {
					Molecule3.push_back(Final[i].xcoord / distconst);
					Molecule3.push_back(Final[i].ycoord / distconst);
				}
				//вектор силы по координате у стал противоположно направленным
				Vectors[i].coordy = (-1) * Vectors[i].coordy;
				Workpack[i].speedy = (-1) * Workpack[i].speedy;
				//оставшееся после столкновения время молекула движется с противоположным направлением по координате у
				Workpack[i].speedx = Workpack[i].speedx + Vectors[i].coordx * t2 / mas;
				Workpack[i].speedy = Workpack[i].speedy + Vectors[i].coordy * t2 / mas;
				Workpack[i].speedz = Workpack[i].speedz + Vectors[i].coordz * t2 / mas;

				Workpack[i].xcoord = Workpack[i].xcoord + Workpack[i].speedx * t2 + Vectors[i].coordx * t2 * t2 / (2 * mas);
				Workpack[i].ycoord = Workpack[i].ycoord + Workpack[i].speedy * t2 + Vectors[i].coordy * t2 * t2 / (2 * mas);
				Workpack[i].zcoord = Workpack[i].zcoord + Workpack[i].speedz * t2 + Vectors[i].coordz * t2 * t2 / (2 * mas);
				Check[i] = 1;
			}
			if (i== 0) {
				Molecule1.push_back(Workpack[i].xcoord / distconst);
			    Molecule1.push_back(Workpack[i].ycoord / distconst);
			}
			if (i == 2) {
				Molecule2.push_back(Workpack[i].xcoord / distconst);
				Molecule2.push_back(Workpack[i].ycoord / distconst);
			}
			if (i == 5) {
				Molecule3.push_back(Workpack[i].xcoord / distconst);
				Molecule3.push_back(Workpack[i].ycoord / distconst);
			}
			
		}
	}
	cout << endl;
	cout << endl;
	cout << endl;
	cout << "Coordinates of collision of molecules with the surface(in nanometers)";
	cout << endl;
	//выводим на экран координаты столкновения молекул с поверхностью, одна единица в координатах соответствует 0.15 нанометрам
	for (int i = 0; i < 20; i++)
	{
		if (Final[i].xcoord != 0)
		{
			cout << "    " << i << "    " << "initial coordinates" << "    " << Initialx[i] / distconst << "    " << Initialy[i] / distconst << "    " << Initialz[i] / distconst << "    " << "final coordinates" << "    " << Final[i].xcoord / distconst << "    " << Final[i].ycoord / distconst << "    " << Final[i].zcoord / distconst << endl;
			cout << endl;
		}
	}
	cout << "Coordinates of molecules which have not collided with the surface after 1 nanosecond" << endl;
	for (int i = 0; i < 20; i++)
	{
		if (Final[i].xcoord == 0)
		{
			cout << endl;

			cout << "    " << i << "    " << "initial coordinates" << "    " << Initialx[i] / distconst << "    " << Initialy[i] / distconst << "    " << Initialz[i] / distconst << "    " << "final coordinates" << "    " << Workpack[i].xcoord / distconst << "    " << Workpack[i].ycoord / distconst << "    " << Workpack[i].zcoord / distconst << endl;
		}
	}
	
	Graphics(Molecule1, Molecule2, Molecule3,1000);
	
	return 0;
}
