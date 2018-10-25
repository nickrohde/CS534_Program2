
#pragma region Includes:
	
	#include <iostream>	// cout, cerr
	#include <chrono>	// timer
	#include <cmath>	// pow

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
	const int DEFAULT_TIME = 500;		// default # of time steps
	const int DEFAULT_INTERVAL = 0;		// default output interval 

	const double WAVE_SPEED = 1.0;		// speed of waves (c)
	const double TIME_QUANTUM = 0.1;	// time steps (dt)
	const double SYS_CHANGE = 2.0;		// change in system (dd)

#pragma endregion


/// <summary>
///          Program entry point, verifies input and calls run.
/// </summary>
/// <param name="argc">
///          The argument count received from operating system.
/// </param>
/// <param name="argv">
///          The argument values received from operating system.
/// </param>
/// <returns>
///          0 on successful execution, 1 otherwise.
/// </returns>
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
	if (!argsValid(size, max_time, interval))
	{
		exit(EXIT_FAILURE);
	} // end if

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


/// <summary>
///          Verifies the values of the command line arguments.
/// </summary>
/// <param name="SIZE">
///          The dimension of the internal square matrices, must be in range [100, 5000].
/// </param>
/// <param name="TIME">
///          The number of time-quanta to simulate, must be in range [3, 5000].
/// </param>
/// <param name="INTERVAL">
///          The interval at which to print out the current wave, must be in range [0, <paramref name="TIME"/>].
/// </param>
/// <returns>
///          True iff SIZE in [100, 5000] AND TIME in [3, 5000] AND INTERVAL in [0, <paramref name="TIME"/>]; otherwise, false.
/// </returns>
bool argsValid(const int SIZE, const int TIME, const int INTERVAL) noexcept
{
	if (!(SIZE >= 100 && SIZE <= 5000 && TIME >= 3 && TIME <= 5000 && INTERVAL >= 0 && INTERVAL <= TIME))
	{
		std::cerr << "usage: Wave2D size max_time interval" << std::endl;
		std::cerr << "       where 100 <= size <= 5000 && 3 <= time <= 5000 && 0 <= interval <= time" << std::endl;
		return false;
	} // end if

	return true;
} // end method argsValid


/// <summary>
///          Allocates the memory of the simulation space.
/// </summary>
/// <param name="z">
///          A pointer to the memory being allocated.
/// </param>
/// <param name="SIZE">
///          The dimension of the internal square matrices.
/// </param>
void instantiate(double***& z, const int SIZE) noexcept
{
	try
	{
		for (int p = 0; p < 3; p++)
		{
			z[p] = new double*[SIZE];
			for (int i = 0; i < SIZE; i++)
			{
				z[p][i] = new double[SIZE];
				memset(z[p][i], 0, SIZE * sizeof(double));
			} // end for i
		} // end for p
	} // end try
	catch (std::bad_alloc& e)
	{
		std::cout << "Memory allocation failed!\n\tReason: " << e.what() << std::endl;
		deallocate(z, SIZE);
		exit(EXIT_FAILURE);
	} // end catch
} // end method instantiate


/// <summary>
///          Deallocates the memory of the simulation space.
/// </summary>
/// <param name="z">
///          A pointer to the memory being deallocated.
/// </param>
/// <param name="SIZE">
///          The dimension of the internal square matrices.
/// </param>
/// <remarks>
///			The value of <paramref name="z"/> will be nullptr after execution.
/// </remarks>
void deallocate(double ***& z, const int SIZE) noexcept
{
	for (int p = 0; p < 3; p++)
	{
		for (int i = 0; i < SIZE; i++)
		{
			delete[] z[p][i];
		} // end for i

		delete[] z[p];
	} // end for p

	delete[] z;
	z = nullptr;
} // end method deallocate


/// <summary>
///          Deallocates the memory of the simulation space.
/// </summary>
/// <param name="z">
///          A pointer to the memory being deallocated.
/// </param>
/// <param name="SIZE">
///          The dimension of the internal square matrices.
/// </param>
/// <remarks>
///			The value of <paramref name="z"/> will be nullptr after execution.
/// </remarks>
void initialize(double** t_0, int weight, const int SIZE) noexcept
{
	// memset has already zeroed out the memory, only need to initialize the wave
	for (auto i = (40 * weight) + 1; i < (60 * weight); i++)
	{
		for (auto j = (40 * weight) + 1; j < (60 * weight); j++)
		{
			t_0[i][j] = 20.0;
		} // end for j
	} // end for i
} // end method initialize


/// <summary>
///          Calculates the wave's state at t = 1.
/// </summary>
/// <param name="t_0">
///          Output parameter, the new state of the wave will be stored in this matrix.
/// </param>
/// <param name="t_1">
///          The state of the wave at t = -1 relative to right now.
/// </param>
/// <param name="SIZE">
///          The dimension of the internal square matrices.
/// </param>
void calculate(double** t_0, double** t_1, const int SIZE) noexcept
{
	for (auto i = 0; i < SIZE; i++)
	{
		for (auto j = 0; j < SIZE; j++)
		{
			if (!i || !j || i + 1 == SIZE || j + 1 == SIZE)
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


/// <summary>
///          Calculates the wave's state at t >= 2.
/// </summary>
/// <param name="t_0">
///          Output parameter, the new state of the wave will be stored in this matrix.
/// </param>
/// <param name="t_1">
///          The state of the wave at t = -1 relative to right now.
/// </param>
/// <param name="t_2">
///          The state of the wave at t = -2 relative to right now.
/// </param>
/// <param name="SIZE">
///          The dimension of the internal square matrices.
/// </param>
void calculate(double** t_0, double** t_1, double** t_2, const int SIZE) noexcept
{
	for (auto i = 0; i < SIZE; i++)
	{
		for (auto j = 0; j < SIZE; j++)
		{
			if (!i || !j || i + 1 == SIZE || j + 1 == SIZE)
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
/// <summary>
///          Outputs the given state of the wave according to the prescribed format.
/// </summary>
/// <param name="current">
///          Current state of the wave for output.
/// </param>
/// <param name="TIME">
///          The time-slice corresponding to this wave state.
/// </param>
/// <param name="SIZE">
///          The dimension of the matrix.
/// </param>
void outputWave(double** current, const int TIME, const int SIZE) noexcept
{
	std::cout << TIME << std::endl;

	for (auto i = 0; i < SIZE; i++)
	{
		for (auto j = 0; j < SIZE; j++)
		{
			std::cout << current[i][j];
		} // end for j

		std::cout << std::endl;
	} // end for i

	std::cout << std::endl;
} // end method outputWave

/*/


//!_START_REMOVE!
/// <summary>
///          DEBUG VERSION 
/// </summary>
void outputWave(double** current, const int TIME, const int SIZE) noexcept
{
	for (auto i = 0; i < SIZE; i++)
	{
		for (auto j = 0; j < SIZE; j++)
		{
			std::cout << current[i][j];
			if (j + 1 < SIZE)
			{
				std::cout << ",";
			}
			else if (i+1 < SIZE)
			{
				std::cout << "|";
			}
		} // end for j
	} // end for i

	std::cout << std::endl;
} // end method outputWave
//!_END_REMOVE!
//*/


/// <summary>
///          Calculates the wave's state at t >= 2.
/// </summary>
/// <param name="SIZE">
///          The dimension of the internal square matrices.
/// </param>
/// <param name="MAX_TIME">
///          The number of time steps to simulate.
/// </param>
/// <param name="INTERVAL">
///          The printing interval for output.
/// </param>
void run(const int SIZE, const int MAX_TIME, const int INTERVAL)
{
	// Variables:
	double*** z = new double**[3], // simulation space
		  **  current = nullptr,
		  **  previous = nullptr,
		  **  previous2 = nullptr,
		  **  temp = nullptr;


	// instantiate memory of z
	instantiate(z, SIZE);

	_TimePoint end, start = _Clock::now();

	// time = 0;
	// initialize the simulation space: calculate z[0][][]
	initialize(z[0], static_cast<int>(SIZE / DEFAULT_SIZE), SIZE);

	// time = 1
	// calculate z[1][][] 
	// cells not on edge
	calculate(z[1], z[0], SIZE);

	current = z[2];
	previous = z[1];
	previous2 = z[0];

	#if !BENCHMARK
		// output starting point
		if (INTERVAL != 0)
		{
			outputWave(z[0], 0, SIZE);
		} // end if
		// output second iteration if needed
		if (INTERVAL == 1)
		{
			outputWave(z[1], 0, SIZE);
		} // end if
	#endif	

	// simulate wave diffusion from time = 2
	for (int t = 2; t < MAX_TIME; t++)
	{
		// calculate new values of wave
		calculate(current, previous, previous2, SIZE);

		#if !BENCHMARK 
			// print out current status every interval iterations
			if (INTERVAL != 0 && !(t % INTERVAL))
			{
				outputWave(current, t, SIZE);
			} // end if
		#endif	

		// rotate pointers
		temp = current;
		current = previous2;
		previous2 = previous;
		previous = temp;
	} // end of simulation

	end = _Clock::now();

	#if !BENCHMARK
		// output end result
		if (INTERVAL != 0)
		{
			outputWave(current, MAX_TIME, SIZE);
		} // end if
	#endif

	std::cout << "Elapsed time: " << DURATION_IN_MICROS(end - start) << " microseconds" << std::endl;

	// null ptrs to memory being freed
	current = previous = previous2 = temp = nullptr;

	deallocate(z, SIZE);
} // end method run


