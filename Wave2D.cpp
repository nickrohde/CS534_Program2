
#pragma region Includes:
	
	#include <iostream>
	#include <chrono>	// timer

#pragma endregion


#pragma region Prototypes:

	// input verification
	bool argsValid(const int size, const int time, const int interval) noexcept;

	// memory management
	void instantiate(double***& z, const int size) noexcept;
	void deallocate(double ***& z, const int size) noexcept;

	void run(const int size, const int max_time, const int interval);

	// simulation methods
	void initialize(double** z, int weight, const int size) noexcept;
	void calculate(double** t_0, double** t_1, const int size) noexcept;
	void calculate(double** t_0, double** t_1, double** t_2, const int size) noexcept;
	void outputWave(double** current, const int t, const int size) noexcept;

#pragma endregion


#pragma region Defines:

	#define DEBUG 0
	#define BENCHMARK 0
	#define DURATION_IN_MICROS(A) std::chrono::duration_cast<std::chrono::microseconds>((A)).count()

#pragma endregion


#pragma region Typedefs:

	typedef std::chrono::high_resolution_clock	_Clock;		// simplify clock name
	typedef _Clock::time_point					_TimePoint;	// simplify time point name

#pragma endregion


#pragma region Global Constants:

	const int DEFAULT_SIZE = 100;		// default system size
	const int DEFAULT_CELL_WIDTH = 8;	// default cell width
	const int DEFAULT_TIME = 20;
	const int DEFAULT_INTERVAL = 1;

	const double WAVE_SPEED = 1.0;		// speed of waves (c)
	const double TIME_QUANTUM = 0.1;	// time steps (dt)
	const double SYS_CHANGE = 2.0;		// change in system (dd)

#pragma endregion


int main(int argc, char **argv)
{
	// Variables:
	int size;
	int max_time;
	int interval;


	// verify arguments
	if (argc == 4) 
	{
		size = atoi(argv[1]);
		max_time = atoi(argv[2]);
		interval = atoi(argv[3]);		
	} // end if
	else if (argc == 3)
	{
		size = atoi(argv[1]);
		max_time = atoi(argv[2]);
		interval = DEFAULT_INTERVAL;
	} // end elif
	else if (argc == 2)
	{
		size = atoi(argv[1]);
		max_time = DEFAULT_TIME;
		interval = DEFAULT_INTERVAL;
	} // end elif
	else
	{
		size = DEFAULT_SIZE;
		max_time = DEFAULT_TIME;
		interval = DEFAULT_INTERVAL;
	} // end else

	// check if args in range
	argsValid(size, max_time, interval) ? __noop : exit(EXIT_FAILURE);

	try
	{
		run(size, max_time, interval);
	} // end try
	catch (std::exception e)
	{
		std::cerr << "An error occurred!\n\tReason: " << e.what() << std::endl;
	} // end catch

	#if (defined(_WIN32) || defined(_WIN64)) && DEBUG
		system("pause");
	#endif

	return 0;
} // end Main


bool argsValid(const int size, const int time, const int interval) noexcept
{
	if (!(size >= 100 && size <= 1000 && time >= 3 && time <= 1000 && interval >= 0 && interval <= 100))
	{
		std::cerr << "usage: Wave2D size max_time interval" << std::endl;
		std::cerr << "       where size >= 100 && time >= 3 && interval >= 0" << std::endl;
		return false;
	} // end if

	return true;
} // end method argsValid


void instantiate(double***& z, const int size) noexcept
{
	try
	{
		for (int p = 0; p < 3; p++)
		{
			z[p] = new double*[size];
			for (int i = 0; i < size; i++)
			{
				z[p][i] = new double[size];
				memset(z[p][i], 0, size * sizeof(double));
			} // end for i
		} // end for p
	} // end try
	catch (std::bad_alloc& e)
	{
		std::cout << "Memory allocation failed!\n\tReason: " << e.what() << std::endl;
		deallocate(z, size);
		exit(EXIT_FAILURE);
	} // end catch
} // end method instantiate


void deallocate(double ***& z, const int size) noexcept
{
	for (int p = 0; p < 3; p++)
	{
		for (int i = 0; i < size; i++)
		{
			delete[] z[p][i];
		} // end for i
		delete[] z[p];
	} // end for p

	delete[] z;
} // end method deallocate


void initialize(double** t_0, int weight, const int size) noexcept
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i > 40 * weight && i < 60 * weight  && j > 40 * weight && j < 60 * weight)
			{
				t_0[i][j] = 20.0;
			} // end if
			else
			{
				t_0[i][j] = 0.0;
			} // end else
		} // end for j
	} // end for i
} // end method initialize


void calculate(double** t_0, double** t_1, const int size) noexcept
{
	for (auto i = 0; i < size; i++)
	{
		for (auto j = 0; j < size; j++)
		{
			if (!i || !j || i + 1 == size || j + 1 == size)
			{
				t_0[i][j] = 0.0;
			} // end if
			else
			{
				t_0[i][j] = t_1[i][j] + (WAVE_SPEED * WAVE_SPEED) / 2.0 * pow((TIME_QUANTUM / SYS_CHANGE), 2.0) * (t_1[i + 1][j] + t_1[i - 1][j] + t_1[i][j + 1] + t_1[i][j - 1] - 4.0 * t_1[i][j]);
			} // end else
		} // end for j
	} // end for i
} // end method calculate(3)


void calculate(double** t_0, double** t_1, double** t_2, const int size) noexcept
{
	for (auto i = 0; i < size; i++)
	{
		for (auto j = 0; j < size; j++)
		{
			if (!i || !j || i + 1 == size || j + 1 == size)
			{
				t_0[i][j] = 0.0;
			} // end if
			else
			{
				t_0[i][j] = 2.0 * t_1[i][j] - t_2[i][j] + (WAVE_SPEED * WAVE_SPEED) * pow((TIME_QUANTUM / SYS_CHANGE), 2.0) * (t_1[i + 1][j] + t_1[i - 1][j] + t_1[i][j + 1] + t_1[i][j - 1] - 4.0 * t_1[i][j]);
			} // end else
		} // end for j
	} // end for i
} // end method calculate(4)


/*

void outputWave(double** current, const int t, const int size) noexcept
{
	std::cout << t << std::endl;

	for (auto i = 0; i < size; i++)
	{
		for (auto j = 0; j < size; j++)
		{
			std::cout << current[i][j];
		} // end for j

		std::cout << std::endl;
	} // end for i

	std::cout << std::endl;
} // end method outputWave

/*/

void outputWave(double** current, const int t, const int size) noexcept
{
	for (auto i = 0; i < size; i++)
	{
		for (auto j = 0; j < size; j++)
		{
			std::cout << current[i][j];
			if (j + 1 < size)
			{
				std::cout << ",";
			}
			else if (i+1 < size)
			{
				std::cout << "|";
			}
		} // end for j
	} // end for i

	std::cout << std::endl;
}
//*/

void run(const int size, const int max_time, const int interval)
{
	// Variables:
	double*** z = new double**[3]; // simulation space
	double** current = nullptr,
		  ** previous = nullptr,
		  ** previous2 = nullptr,
		  ** temp = nullptr;

	// instantiate memory of z
	instantiate(z, size);

	_TimePoint end, start = _Clock::now();

	// time = 0;
	// initialize the simulation space: calculate z[0][][]
	initialize(z[0], static_cast<int>(size / DEFAULT_SIZE), size);

	// time = 1
	// calculate z[1][][] 
	// cells not on edge
	calculate(z[1], z[0], size);

	current = z[2];
	previous = z[1];
	previous2 = z[0];

	#if !BENCHMARK
		// output starting point
		outputWave(z[0], 0, size);
		
		// output second iteration if needed
		if (interval <= 1)
		{
			outputWave(z[1], 0, size);
		} // end if
	#endif	

	// simulate wave diffusion from time = 2
	for (int t = 2; t < max_time; t++)
	{
		// calculate new values of wave
		calculate(current, previous, previous2, size);

		#if !BENCHMARK 
		// print out current status every interval iterations
			if (t <= 1 || !(t % interval))
			{
				outputWave(current, t, size);
			} // end if
		#endif	

		// rotate pointers
		temp = current;
		current = previous2;
		previous2 = previous;
		previous = temp;
	} // end of simulation

	end = _Clock::now();

	std::cout << "Elapsed time: " << DURATION_IN_MICROS(end - start) << std::endl;

	// null ptrs to memory being freed
	current = previous = previous2 = temp = nullptr;

	deallocate(z, size);
} // end method run


