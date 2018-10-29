
#pragma region Includes:
	
	#include <iostream>		// cout, cerr
	#include <chrono>		// timer
	#include <cmath>		// pow
	#include <mpi.h>		// MPI functionality
	#include <vector>		// vector
	#include <thread>		// debug
	#include <sstream>		// debug
	#include <fstream>		// ofstream
	#include <omp.h>		// OpenMP directives
	#include <string.h>		// memset

#pragma endregion


#pragma region Defines:

	#define DEBUG 0
	#define DEBUG2 1
	#define BENCHMARK 0
	#define __WORLD MPI_COMM_WORLD
	#define DURATION_IN_MICROS(A) std::chrono::duration_cast<std::chrono::microseconds>((A)).count()

	#ifndef EXIT_FAILURE
		#define EXIT_FAILURE 1
	#endif
	#ifndef EXIT_SUCCESS
		#define EXIT_SUCCESS 0
	#endif

#pragma endregion


#pragma region Typedefs:

	typedef std::chrono::high_resolution_clock	_Clock;		// simplify clock name
	typedef _Clock::time_point					_TimePoint;	// simplify time point name

#pragma endregion


#pragma region Structures:
	
	struct Comm_Info
	{
		double* left_recv;
		double* left_send;
		double* right_recv;
		double* right_send;
		int size;
	}; // end struct Comm_Info

#pragma endregion


#pragma region Prototypes:

	// DEBUG
	void run2(const int SIZE, const int MAX_TIME, const int INTERVAL);
	void calculate2(double** t_0, double** t_1, const int COLS, const int ROWS) ;
	void dumpMatrix(double** t_0, const int COLS, const int ROWS);

	// input verification
	bool argsValid(const int SIZE, const int TIME, const int INTERVAL, const int THREADS);
	void processArgs(const char** argv, const int argc, int* size, int* max_time, int* interval, int* threads);

	// memory management
	void instantiate(double***& z, const int COLS, const int ROWS) ;
	void deallocate(double ***& z, const int COLS) ;

	// MPI functions
	//std::vector<MPI_Request*>* sendInitial(double** current, const int SIZE, const int MPI_SIZE);
	void sendInitial(double** current, const int SIZE);
	void waitAll(std::vector<MPI_Request*>*& requests);
	void receiveInitial(double** current, const int COLS, const int ROWS);
	void setSyncInfo(Comm_Info& comm_info, double** current, const int COLS, const int ROWS);
	void communicateOverlap(Comm_Info& comm_info);
	void reconstructMatrix(double** current, const int COLS, const int ROWS);

	void run(const int SIZE, const int TIME, const int INTERVAL);

	// simulation methods
	void initialize(double** z, int weight, const int COLS, const int ROWS) ;
	void calculate(double** t_0, double** t_1, const int COLS, const int ROWS) ;
	void calculate(double** t_0, double** t_1, double** t_2, const int COLS, const int ROWS) ;

	void outputWave(double** current, const int TIME, const int COLS, const int ROWS) ;

#pragma endregion


#pragma region Global Constants:

	const int DEFAULT_SIZE = 100;		// default system size
	const int DEFAULT_TIME = 500;		// default # of time steps
	const int DEFAULT_INTERVAL = 0;		// default output interval
	const int DEFAULT_THREADS = 1;

	const double WAVE_SPEED = 1.0;		// speed of waves (c)
	const double TIME_QUANTUM = 0.1;	// time steps (dt)
	const double SYS_CHANGE = 2.0;		// change in system (dd)

#pragma endregion


#pragma region Globals:

	int my_rank;
	int mpi_size;
	int threads;

#pragma endregion


#pragma region Main Program:

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
		int size,
			max_time,
			interval;

		
		try
		{
			// process the command line arguments received from OS
			processArgs(argv, argc, &size, &max_time, &interval, &threads);

			// initialize MPI and OMP
			MPI_Init(&argc, &argv);
			MPI_Bcast(&threads, 1, MPI_INT, 0, __WORLD);
			omp_set_num_threads(threads);

			// run the Wave simmulation
			run(size, max_time, interval);
		} // end try
		catch (std::exception e)
		{
			std::cerr << "An error occurred!\n\tReason: " << e.what() << std::endl;
			exit(EXIT_FAILURE);
		} // end catch

		#if (defined(_WIN32) || defined(_WIN64)) && DEBUG
			system("pause");
		#endif

		return EXIT_SUCCESS;
	} // end Main


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

		int cols,
			rows;

		// start the clock for the simulation
		_TimePoint end, 
			       start = _Clock::now();

		Comm_Info sync_data;


		// get rank and number of nodes
		MPI_Comm_rank(__WORLD, &my_rank);
		MPI_Comm_size(__WORLD, &mpi_size);

		// calculate slice size
		cols = (SIZE / mpi_size) + ((SIZE % mpi_size > my_rank) ? 1 : 0);
		rows = SIZE;

		if (!my_rank)
		{
			start = _Clock::now();
			// instantiate memory of z
			instantiate(z, SIZE + 2, rows);
			// initialize the simulation space: calculate z[0][][]
			initialize(z[0], static_cast<int>(SIZE / DEFAULT_SIZE), SIZE, SIZE);
			// send out initial data to all ranks
			sendInitial(z[0], SIZE);
		} // end if
		else
		{
			// allocate the memory
			instantiate(z, cols + 2, rows);
			// receive slice of work from the manager rank
			receiveInitial(z[0], cols, rows);
		} // end else

		// calculate t = 1
		calculate(z[1], z[0], cols, rows);

		// communicate the boundary data to neighbors
		setSyncInfo(sync_data, z[1], cols, rows);
		communicateOverlap(sync_data);

		// initialize pointers
		current = z[2];
		previous = z[1];
		previous2 = z[0];

	#if !BENCHMARK // ommit this output code in benchmark mode
		// output starting point
		if (INTERVAL != 0)
		{
			if (!my_rank)
			{
				outputWave(z[0], 0, SIZE, SIZE);
			} // end if

			  // output second iteration if needed
			if (INTERVAL == 1)
			{
				if (!my_rank)
				{
					reconstructMatrix(previous, SIZE, SIZE);
					outputWave(z[1], 0, SIZE, SIZE);
				} // end if
				else
				{
					reconstructMatrix(previous, cols, rows);
				} // end else
			} // end if
		} // end if

	#endif	

		// simulate wave diffusion from time = 2
		for (int t = 2; t < MAX_TIME; t++)
		{
			// calculate new values of wave
			calculate(current, previous, previous2, cols, rows);

			// prepare for communication
			setSyncInfo(sync_data, current, cols, rows);

			// send boundary data to neighbours
			communicateOverlap(sync_data);

			#if !BENCHMARK 
				// print out current status
				if (INTERVAL != 0 && !(t % INTERVAL))
				{
					if (!my_rank)
					{
						reconstructMatrix(current, SIZE, SIZE);
						outputWave(current, t, SIZE, SIZE);
					} // end if
					else
					{
						reconstructMatrix(current, cols, rows);
					} // end else
				} // end if
			#endif

			// rotate pointers
			temp = current;
			current = previous2;
			previous2 = previous;
			previous = temp;
		} // end of simulation

		// stop the clock
		end = _Clock::now();

	#if !BENCHMARK
		// output end result
		if (INTERVAL > 1)
		{
			if (!my_rank)
			{
				reconstructMatrix(previous, SIZE, SIZE);
				outputWave(previous, MAX_TIME, SIZE, SIZE);
			} // end if
			else
			{
				reconstructMatrix(previous, cols, rows);
			} // end else
		} // end if
	#endif

		// let MPI know we're done
		MPI_Finalize(); // moving this statement causes MPI to not shut down properly

		// output final result to file
		if (!my_rank)
		{
			std::ofstream results("run_info.txt", std::ios::app | std::ios::out);

			if (!results.is_open() || results.bad())
			{
				std::cerr << "Opening of results file failed! Results will not be stored!\n";
			} // end if
			else
			{
				int start_index = 0,
					end_index = (SIZE / mpi_size) + ((SIZE % mpi_size > 0) ? 1 : 0);

				results << "Execution with " << mpi_size << " MPI nodes and " << threads << " threads.\n";

				for (auto i = 0; i < mpi_size; i++)
				{
					results << "rank[" << i << "]'s range = " << start_index << " ~ " << end_index << std::endl;
					start_index += (SIZE / mpi_size) + ((SIZE % mpi_size > i) ? 1 : 0);
					end_index += (SIZE / mpi_size) + ((SIZE % mpi_size > i) ? 1 : 0);
				} // end for

				results << "Elapsed time: " << DURATION_IN_MICROS(end - start) << " microseconds\n\n";
			} // end else

			results.close();
		} // end if

		// destroy pointers to memory being freed
		current = previous = previous2 = temp = nullptr;

		// free memory used for the wave
		deallocate(z, SIZE);
	} // end method run
	
		
	/// <summary>
	///          Processes the arguments from argv and verifies the values received.
	/// </summary>
	/// <param name="argv">
	///          The command line argument values.
	/// </param>
	/// <param name="argc">
	///          The number of arguments stored in argv.
	/// </param>
	/// <param name="size">
	///          Pointer to the variable where the size should be stored.
	/// </param>
	/// <param name="max_time">
	///          Pointer to the variable where the number of iterations should be stored.
	/// </param>
	/// <param name="interval">
	///          Pointer to the variable where the print interval should be stored.
	/// </param>
	/// <param name="threads">
	///          Pointer to the variable where the # of threads should be stored.
	/// </param>
	void processArgs(const char** argv, const int argc, int* size, int* max_time, int* interval, int* threads)
	{
		// set all arguments to default in case any were ommitted
		*size = DEFAULT_SIZE,
		*max_time = DEFAULT_TIME,
		*interval = DEFAULT_INTERVAL;
		*threads = DEFAULT_THREADS;

		// process command line args
		switch (argc)
		{
		case 5:
			*threads = atoi(argv[4]);

		case 4:
			*interval = atoi(argv[3]);

		case 3:
			*max_time = atoi(argv[2]);

		case 2:
			*size = atoi(argv[1]);
			break;

		default:
			std::cerr << "No arguments received, using defaults." << std::endl;
			break;
		} // end switch

		// verify that the arguments received are actually valid
		if (!argsValid(*size, *max_time, *interval, *threads))
		{
			exit(EXIT_FAILURE);
		} // end if
	} // end method processArgs


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
	///          True iff SIZE >= 100 AND TIME >= 3 AND INTERVAL in >= 0 AND THREADS >= 0; otherwise, false.
	/// </returns>
	bool argsValid(const int SIZE, const int TIME, const int INTERVAL, const int THREADS)
	{
		if (!(SIZE >= 100 && TIME >= 3 && INTERVAL >= 0 && THREADS >= 1))
		{
			std::cerr << "usage: Wave2D size max_time interval num_threads" << std::endl;
			std::cerr << "       where 100 <= size && 3 <= time && 0 <= interval && 1 <= num_threads" << std::endl;
			return false;
		} // end if

		return true;
	} // end method argsValid

#pragma endregion


#pragma region Memory Management:

	/// <summary>
	///          Allocates the memory of the simulation space.
	/// </summary>
	/// <param name="z">
	///          A pointer to the memory being allocated.
	/// </param>
	/// <param name="COLS">
	///          The outer dimension of the matrix.
	/// </param>
	/// <param name="ROWS">
	///          The inner dimension of the matrix.
	/// </param>
	void instantiate(double***& z, const int COLS, const int ROWS) 
	{
		try
		{
			#pragma omp parallel for
			for (int p = 0; p < 3; p++)
			{
				z[p] = new double*[COLS];
				for (int i = 0; i < COLS; i++)
				{
					z[p][i] = new double[ROWS];
					memset(z[p][i], 0, ROWS * sizeof(double));
				} // end for i
			} // end for p
		} // end try
		catch (std::bad_alloc& e)
		{
			std::cout << "Memory allocation failed!\n\tReason: " << e.what() << std::endl;
			deallocate(z, COLS);
			exit(EXIT_FAILURE);
		} // end catch
	} // end method instantiate


	/// <summary>
	///          Deallocates the memory of the simulation space.
	/// </summary>
	/// <param name="z">
	///          A pointer to the memory being deallocated.
	/// </param>
	/// <param name="COLS">
	///          The outer dimension of the matrix.
	/// </param>
	void deallocate(double ***& z, const int COLS) 
	{
		for (int p = 0; p < 3; p++)
		{
			for (int i = 0; i < COLS; i++)
			{
				delete[] z[p][i];
			} // end for i

			delete[] z[p];
		} // end for p

		delete[] z;
		z = nullptr;
	} // end method deallocate

#pragma endregion


#pragma region Data Management:

	/// <summary>
	///          Initializes the simulation space by placing a central tidal wave.
	/// </summary>
	/// <param name="t_0">
	///          Output parameter where the initial wave will be stored.
	/// </param>
	/// <param name="weight">
	///          Constant controlling the size of the initial wave.
	/// </param>
	/// <param name="COLS">
	///          The outer dimension of the matrix.
	/// </param>
	/// <param name="ROWS">
	///          The inner dimension of the matrix.
	/// </param>
	void initialize(double** t_0, int weight, const int COLS, const int ROWS) 
	{
		// memory is already zeroed, just need to place the initial wave
		for (auto i = static_cast<std::size_t>(0.4 * COLS) + 1; (i < (0.6 * COLS)) && (i < COLS); i++)
		{
			for (auto j = static_cast<std::size_t>(0.4 * ROWS) + 1; (j < (0.6 * ROWS)) && (j < ROWS); j++)
			{
				t_0[i][j] = 20.0;
			} // end for j
		} // end for i
	} // end method initialize


	/// <summary>
	///          Outputs the given state of the wave according to the prescribed format.
	/// </summary>
	/// <param name="current">
	///          Current state of the wave for output.
	/// </param>
	/// <param name="TIME">
	///          The time-slice corresponding to this wave state.
	/// </param>
	/// <param name="COLS">
	///          The outer dimension of the matrix.
	/// </param>
	/// <param name="ROWS">
	///          The inner dimension of the matrix.
	/// </param>
	void outputWave(double** current, const int TIME, const int COLS, const int ROWS)
	{
		std::cout << TIME << std::endl;

		for (auto i = 1; i <= COLS; i++)
		{
			for (auto j = 0; j < ROWS; j++)
			{
				std::cout << current[i][j] << " ";
			} // end for j

			std::cout << std::endl;
		} // end for i

		std::cout << std::endl;
	} // end method outputWave

#pragma endregion


/// <summary>
///          Calculates the wave's state at t = 1.
/// </summary>
/// <param name="t_0">
///          Output parameter, the new state of the wave will be stored in this matrix.
/// </param>
/// <param name="t_1">
///          The state of the wave at t = -1 relative to right now.
/// </param>
/// <param name="COLS">
///          The outer dimension of the matrix.
/// </param>
/// <param name="ROWS">
///          The inner dimension of the matrix.
/// </param>
void calculate(double** t_0, double** t_1, const int COLS, const int ROWS) 
{
	#pragma omp parallel for schedule(guided)
	for (auto i = 1; i <= COLS; i++)
	{
		for (auto j = 0; j < ROWS; j++)
		{
			// manager has left boundary
			if (!my_rank)
			{
				if (i == 1 || !j || j + 1 == ROWS)
				{
					t_0[i][j] = 0.0;
				} // end if
				else
				{
					t_0[i][j] = t_1[i][j] + (WAVE_SPEED * WAVE_SPEED) / 2.0 * pow((TIME_QUANTUM / SYS_CHANGE), 2.0) * (t_1[i + 1][j] + t_1[i - 1][j] + t_1[i][j + 1] + t_1[i][j - 1] - 4.0 * t_1[i][j]);
				} // end else
			} // end if
			else if (my_rank + 1 == mpi_size) // right-most worker has right boundary
			{
				if (!j || i + 2 == COLS || j + 1 == ROWS)
				{
					t_0[i][j] = 0.0;
				} // end if
				else
				{
					t_0[i][j] = t_1[i][j] + (WAVE_SPEED * WAVE_SPEED) / 2.0 * pow((TIME_QUANTUM / SYS_CHANGE), 2.0) * (t_1[i + 1][j] + t_1[i - 1][j] + t_1[i][j + 1] + t_1[i][j - 1] - 4.0 * t_1[i][j]);
				} // end else
			} // end elif
			else // all other ranks have no right and left boundaries
			{
				if (!j || j + 1 == ROWS)
				{
					t_0[i][j] = 0.0;
				} // end if
				else
				{
					t_0[i][j] = t_1[i][j] + (WAVE_SPEED * WAVE_SPEED) / 2.0 * pow((TIME_QUANTUM / SYS_CHANGE), 2.0) * (t_1[i + 1][j] + t_1[i - 1][j] + t_1[i][j + 1] + t_1[i][j - 1] - 4.0 * t_1[i][j]);
				} // end else
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
/// <param name="COLS">
///          The outer dimension of the matrix.
/// </param>
/// <param name="ROWS">
///          The inner dimension of the matrix.
/// </param>
void calculate(double** t_0, double** t_1, double** t_2, const int COLS, const int ROWS) 
{
	#pragma omp parallel for schedule(guided)
	for (auto i = 1; i <= COLS; i++)
	{
		for (auto j = 0; j < ROWS; j++)
		{
			// manager has left boundary
			if (!my_rank)
			{
				if (i == 1 || !j || j + 1 == ROWS)
				{
					t_0[i][j] = 0.0;
				} // end if
				else
				{
					t_0[i][j] = 2.0 * t_1[i][j] - t_2[i][j] + (WAVE_SPEED * WAVE_SPEED) * pow((TIME_QUANTUM / SYS_CHANGE), 2.0) * (t_1[i + 1][j] + t_1[i - 1][j] + t_1[i][j + 1] + t_1[i][j - 1] - 4.0 * t_1[i][j]);
				} // end else
			} // end if
			else if (my_rank + 1 == mpi_size) // right-most worker has right boundary
			{
				if (!j || i + 2 == COLS || j + 1 == ROWS)
				{
					t_0[i][j] = 0.0;
				} // end if
				else
				{
					t_0[i][j] = 2.0 * t_1[i][j] - t_2[i][j] + (WAVE_SPEED * WAVE_SPEED) * pow((TIME_QUANTUM / SYS_CHANGE), 2.0) * (t_1[i + 1][j] + t_1[i - 1][j] + t_1[i][j + 1] + t_1[i][j - 1] - 4.0 * t_1[i][j]);
				} // end else
			} // end elif
			else // all other ranks have no right and left boundaries
			{
				if (!j || j + 1 == ROWS)
				{
					t_0[i][j] = 0.0;
				} // end if
				else
				{
					t_0[i][j] = 2.0 * t_1[i][j] - t_2[i][j] + (WAVE_SPEED * WAVE_SPEED) * pow((TIME_QUANTUM / SYS_CHANGE), 2.0) * (t_1[i + 1][j] + t_1[i - 1][j] + t_1[i][j + 1] + t_1[i][j - 1] - 4.0 * t_1[i][j]);
				} // end else
			} // end else			
		} // end for j
	} // end for i
} // end method calculate(4)


#pragma region Communication:

	/// <summary>
	///          Sends the initial data from the manager node to all mpi workers.
	/// </summary>
	/// <param name="current">
	///          The current matrix descriptor for the wave.
	/// </param>
	/// <param name="SIZE">
	///          The dimensions of the square matrix.
	/// </param>
	//std::vector<MPI_Request*>* sendInitial(double** current, const int SIZE, const int MPI_SIZE)
	void sendInitial(double** current, const int SIZE)
	{
		if (mpi_size == 1) { return; }

		const int SLICE_SIZE = SIZE / mpi_size;
		const bool UNEVEN = (SIZE % mpi_size) > 0;

		//auto requests = new std::vector<MPI_Request*>();
		//requests->reserve(SIZE);

		// special case for right boundary worker
		// the right boundary is guaranteed to never get extra work
		if (DEBUG && std::cerr << "rank 0 assigned 0 ~ " << SLICE_SIZE - 1 << std::endl);
		for (auto i = ((mpi_size - 1) * (SLICE_SIZE)); i <= SIZE; i++)
		{
			//std::cout << "sending row " << i << " to " << mpi_size - 1 << std::endl;
			MPI_Send(current[i], SIZE, MPI_DOUBLE, mpi_size - 1, 0, __WORLD);
		} // end for i
	
		// send data to remaining nodes
		for (auto dest = 1; dest < mpi_size - 1; dest++)
		{
			int offset = (SIZE % mpi_size > dest) ? 1 : 0;

			if (DEBUG && std::cerr << "rank " << dest << " assigned " << (SLICE_SIZE * dest) - 1 << " ~ " << SLICE_SIZE * (dest + 1) + 1 << std::endl);
			for (auto i = (SLICE_SIZE * dest); (i <= SIZE) && (i <= SLICE_SIZE * (dest + 1) + 1); i++)
			{
				//std::cout << "sending row " << i << " to " << dest << std::endl;
				MPI_Send(current[i], SIZE, MPI_DOUBLE, dest, 0, __WORLD);
			} // end for i
		} // end for dest

		if (DEBUG && std::cerr << "rank " << (mpi_size - 1) << " assigned " << ((mpi_size - 1) * (SLICE_SIZE)-1) << " ~ " << SIZE << std::endl);

		//return requests;
	} // end method sendInitial


	/// <summary>
	///          Receives the initial data from the manager node on the mpi worker nodes.
	/// </summary>
	/// <param name="current">
	///          The current matrix descriptor for the wave where data will be placed.
	/// </param>
	/// <param name="COLS">
	///          The number of columns this worker is in charge of.
	/// </param>
	/// <param name="ROWS">
	///          The number of items per row.
	/// </param>
	void receiveInitial(double** current, const int COLS, const int ROWS)
	{
		if (mpi_size == 1) { return; }

		if (my_rank == mpi_size - 1)
		{
			// receive matrix slices from manager
			for (auto i = 0; i <= COLS; i++)
			{
				MPI_Recv(current[i], ROWS, MPI_DOUBLE, 0, 0, __WORLD, MPI_STATUS_IGNORE);
			} // end for i
		} // end if
		else
		{
			// receive matrix slices from manager
			for (auto i = 0; i <= COLS + 1; i++)
			{
				MPI_Recv(current[i], ROWS, MPI_DOUBLE, 0, 0, __WORLD, MPI_STATUS_IGNORE);
			} // end for i
		} // end else

	
	} // end method receiveInitial


	/// <summary>
	///          Communicates the overlapping data between neighboring mpi ranks.
	/// </summary>
	/// <param name="comm_info">
	///          Structure containing pointers to the data to send, data to receive, and number of items per row.
	/// </param>
	void communicateOverlap(Comm_Info& comm_info)
	{
		if (mpi_size == 1) { return; }
		//auto requests = new std::vector<MPI_Request*>();
		//std::stringstream ss;

		//ss << "com_log_rank_" << comm_info.my_rank << ".log";

		//std::string name = ss.str();

		//std::ofstream logfile(name, std::ios::out | std::ios::app);

		//logfile << "[" << comm_info.my_rank << "]: entered communicateOverlap." << std::endl;

		if (my_rank == mpi_size - 1) // right slice 
		{
			//if (DEBUG && logfile << "rank (right boundary) " << comm_info.my_rank << " preparing to communicate with " << comm_info.my_rank - 1 << std::endl);
			int parity = my_rank % 2;
			auto sent = false;

			for (auto i = 0; i < 2; i++)
			{
				if ((parity || i) && !sent) // odd rank or has already received data
				{
					//if (DEBUG && logfile << "rank (right boundary) " << comm_info.my_rank << " sending left boundary to " << comm_info.my_rank - 1 << std::endl);
					//auto temp = new MPI_Request();
					//MPI_Isend(comm_info.left_send, comm_info.size, MPI_DOUBLE, comm_info.my_rank - 1, 2, __WORLD, temp);
					MPI_Send(comm_info.left_send, comm_info.size, MPI_DOUBLE, my_rank - 1, 2, __WORLD);
					//requests->push_back(temp);
					sent = true;
				} // end if
				else // even rank or has already sent data
				{
					//if( DEBUG && logfile << "rank (right boundary) " << comm_info.my_rank << " receiving left boundary from " << comm_info.my_rank - 1 << std::endl);
					//auto temp = new MPI_Request();
					//MPI_Irecv(comm_info.left_recv, comm_info.size, MPI_DOUBLE, comm_info.my_rank - 1, 2, __WORLD, temp);
					MPI_Recv(comm_info.left_recv, comm_info.size, MPI_DOUBLE, my_rank - 1, 2, __WORLD, MPI_STATUS_IGNORE);
					//requests->push_back(temp);
				} // end else
			} // end for i
		} // end if
		else if (!my_rank) // left slice (manager)
		{
			//auto temp = new MPI_Request();

			//if (DEBUG && logfile << "rank (manager) " << comm_info.my_rank << " receiving right boundary from " << 1 << std::endl);
			//MPI_Irecv(comm_info.right_recv, comm_info.size, MPI_DOUBLE, 1, 2, __WORLD, temp);
			MPI_Recv(comm_info.right_recv, comm_info.size, MPI_DOUBLE, 1, 2, __WORLD, MPI_STATUS_IGNORE);
			//requests->push_back(temp);

			//temp = new MPI_Request();
			//if (DEBUG && logfile << "rank (manager) " << comm_info.my_rank << " sending right boundary to " << 1 << std::endl);
			//MPI_Isend(comm_info.right_send, comm_info.size, MPI_DOUBLE, 1, 2, __WORLD, temp);
			MPI_Send(comm_info.right_send, comm_info.size, MPI_DOUBLE, 1, 2, __WORLD);
			//requests->push_back(temp);
		} // end elif	
		else if (my_rank % 2) // odd slices
		{
			//auto temp = new MPI_Request();
			//if (DEBUG && logfile << "rank " << comm_info.my_rank << " sending right boundary to " << comm_info.my_rank + 1 << std::endl);
			//MPI_Isend(comm_info.right_send, comm_info.size, MPI_DOUBLE, comm_info.my_rank + 1, 2, __WORLD, temp);
			MPI_Send(comm_info.right_send, comm_info.size, MPI_DOUBLE, my_rank + 1, 2, __WORLD);
			//requests->push_back(temp);

			//temp = new MPI_Request();
			//if (DEBUG && logfile << "rank " << comm_info.my_rank << " sending left boundary to " << comm_info.my_rank - 1 << std::endl);
			//MPI_Isend(comm_info.left_send, comm_info.size, MPI_DOUBLE, comm_info.my_rank - 1, 2, __WORLD, temp);
			MPI_Send(comm_info.left_send, comm_info.size, MPI_DOUBLE, my_rank - 1, 2, __WORLD);
			//requests->push_back(temp);

			//temp = new MPI_Request();
			//if (DEBUG && logfile << "rank " << comm_info.my_rank << " receiving left boundary from " << comm_info.my_rank - 1 << std::endl);
			//MPI_Irecv(comm_info.left_recv, comm_info.size, MPI_DOUBLE, comm_info.my_rank - 1, 2, __WORLD, temp);
			MPI_Recv(comm_info.left_recv, comm_info.size, MPI_DOUBLE, my_rank - 1, 2, __WORLD, MPI_STATUS_IGNORE);
			//requests->push_back(temp);

			//temp = new MPI_Request();
			//if (DEBUG && logfile << "rank " << comm_info.my_rank << " receiving right boundary from " << comm_info.my_rank + 1 << std::endl);
			//MPI_Irecv(comm_info.right_recv, comm_info.size, MPI_DOUBLE, comm_info.my_rank + 1, 2, __WORLD, temp);
			MPI_Recv(comm_info.right_recv, comm_info.size, MPI_DOUBLE, my_rank + 1, 2, __WORLD, MPI_STATUS_IGNORE);
			//requests->push_back(temp);

		} // end elif	
		else // even slices
		{
			//auto temp = new MPI_Request();
			//if (DEBUG && logfile << "rank " << comm_info.my_rank << " receiving left boundary from " << comm_info.my_rank - 1 << std::endl);
			//MPI_Irecv(comm_info.left_recv, comm_info.size, MPI_DOUBLE, comm_info.my_rank - 1, 2, __WORLD, temp);
			MPI_Recv(comm_info.left_recv, comm_info.size, MPI_DOUBLE, my_rank - 1, 2, __WORLD, MPI_STATUS_IGNORE);
			//requests->push_back(temp);

			//temp = new MPI_Request();
			//if (DEBUG && logfile << "rank " << comm_info.my_rank << " receiving right boundary from " << comm_info.my_rank + 1 << std::endl);
			//MPI_Irecv(comm_info.right_recv, comm_info.size, MPI_DOUBLE, comm_info.my_rank + 1, 2, __WORLD, temp);
			MPI_Recv(comm_info.right_recv, comm_info.size, MPI_DOUBLE, my_rank + 1, 2, __WORLD, MPI_STATUS_IGNORE);
			//requests->push_back(temp);
		

			//temp = new MPI_Request();
			//if (DEBUG && logfile << "rank " << comm_info.my_rank << " sending right boundary to " << comm_info.my_rank + 1 << std::endl);
			//MPI_Isend(comm_info.right_send, comm_info.size, MPI_DOUBLE, comm_info.my_rank + 1, 2, __WORLD, temp);
			MPI_Send(comm_info.right_send, comm_info.size, MPI_DOUBLE, my_rank + 1, 2, __WORLD);
			//requests->push_back(temp);

			//temp = new MPI_Request();
			//if (DEBUG && logfile << "rank " << comm_info.my_rank << " sending left boundary to " << comm_info.my_rank - 1 << std::endl);
			//MPI_Isend(comm_info.left_send, comm_info.size, MPI_DOUBLE, comm_info.my_rank - 1, 2, __WORLD, temp);
			MPI_Send(comm_info.left_send, comm_info.size, MPI_DOUBLE, my_rank - 1, 2, __WORLD);
			//requests->push_back(temp);
		} // end else

		//logfile.close();
	} // end method communicateOverlap


	/// <summary>
	///          Reconstructs the full matrix on the manager node for output.
	/// </summary>
	/// <param name="current">
	///          The current matrix descriptor for the wave that will be reconstructed from the slices.
	/// </param>
	/// <param name="COLS">
	///          The number of columns in the matrix.
	/// </param>
	/// <param name="ROWS">
	///          The number of items per row.
	/// </param>
	void reconstructMatrix(double** current, const int COLS, const int ROWS)
	{
		if (mpi_size == 1) { return; }
		//auto requests = new std::vector<MPI_Request*>();
		//requests->reserve(COLS);

		// Manager receives data from all ranks
		if (!my_rank)
		{
			const auto SLICE_SIZE = COLS / mpi_size;

			for (auto src = 1; src < mpi_size; src++)
			{
				for (auto i = (SLICE_SIZE * src) + 1; (i <= COLS) && (i <= SLICE_SIZE * (src + 1)); i++)
				{
					if (DEBUG && std::cerr << "rank " << my_rank << " receiving " << i << " from " << src << std::endl);
					//auto temp = new MPI_Request();
					//MPI_Irecv(current[i], ROWS, MPI_DOUBLE, src, 1, __WORLD, temp);
					MPI_Recv(current[i], ROWS, MPI_DOUBLE, src, 1, __WORLD, MPI_STATUS_IGNORE);
					//requests->push_back(temp);
				} // end for i
			} // end for dest
		} // end if
		else // Workers send data to manager  
		{
			for (auto i = 1; i <= COLS; i++)
			{
				if (DEBUG && std::cerr << "rank " << my_rank << " sending " << i << " to manager" << std::endl);
				//auto temp = new MPI_Request();
				//MPI_Isend(current[i], ROWS, MPI_DOUBLE, 0, 1, __WORLD, temp);
				MPI_Send(current[i], ROWS, MPI_DOUBLE, 0, 1, __WORLD);
				//requests->push_back(temp);
			} // end for i
		} // end else

		//waitAll(requests);
		//requests->clear();
		//delete requests;

		//requests = nullptr;
	} // end method reconstructMatrix


	/// <summary>
	///          Prepares the Comm_Info object for communicating the overlap.
	/// </summary>
	/// <param name="current">
	///          The current matrix descriptor for the wave.
	/// </param>
	/// <param name="COLS">
	///          The number of columns in the matrix.
	/// </param>
	/// <param name="ROWS">
	///          The number of items per row.
	/// </param>
	void setSyncInfo(Comm_Info& comm_info, double** current, const int COLS, const int ROWS)
	{
		comm_info.left_recv = current[0];
		comm_info.right_recv = current[COLS+1];

		comm_info.left_send = current[1];
		comm_info.right_send = current[COLS];

		comm_info.size = ROWS;
	} // end method setSyncInfo


	/// <summary>
	///          Waits for all requests in the <paramref name="requests"/> to be completed and destroys the objects.
	/// </summary>
	/// <param name="requests">
	///          A vector containing MPI requests that should be waited for.
	/// </param>
	void waitAll(std::vector<MPI_Request*>*& requests)
	{
		if (mpi_size == 1) { return; }
		for (auto& req : *requests)
		{
			MPI_Status stat;
			MPI_Wait(req, &stat);

			if (stat.MPI_ERROR != MPI_SUCCESS)
			{
				std::cerr << "MPI Request failed!\n\tCode: " << stat.MPI_ERROR << "\n\tFrom: " << stat.MPI_SOURCE << "\n";
			} // end if
			delete req;
		} // end for
	} // end method synchronize

#pragma endregion


#pragma region Debug Code

#if DEBUG

	void run2(const int SIZE, const int MAX_TIME, const int INTERVAL)
	{
		double*** z = new double**[3];

		int cols,
			rows;

		_TimePoint end, start = _Clock::now();

		//std::vector<MPI_Request*>* requests;
		Comm_Info sync_data;

		// get rank and number of nodes
		MPI_Comm_rank(__WORLD, &my_rank);
		MPI_Comm_size(__WORLD, &mpi_size);

		cols = (SIZE / mpi_size);
		rows = SIZE;

		if (!my_rank)
		{
			// instantiate memory of z
			instantiate(z, SIZE + 2, rows);

			for (auto j = 1; j <= SIZE; j++)
			{
				for (auto k = 0; k < SIZE; k++)
				{
					z[0][j][k] = k + (j - 1) * SIZE;
				}
			}

			// send slices of initial matrix to all processes
			sendInitial(z[0], SIZE);
			dumpMatrix(z[0], cols, rows);
			dumpMatrix(z[0], SIZE, rows);
		} // end if
		else
		{
			instantiate(z, cols + 2, rows);
			receiveInitial(z[0], cols, rows);
			dumpMatrix(z[0], cols, rows);
		} // end else

		if (!my_rank)
		{
			for (auto i = 1; i <= SIZE; i++)
			{
				for (auto j = 0; j < SIZE; j++)
				{
					if (z[0][i][j] < 10)
					{
						std::cout << z[0][i][j] << "    ";
					}
					else if (z[0][i][j] < 100)
					{
						std::cout << z[0][i][j] << "   ";
					}
					else if (z[0][i][j] < 1000)
					{
						std::cout << z[0][i][j] << "  ";
					}
					else
					{
						std::cout << z[0][i][j] << " ";
					}
				}
				std::cout << std::endl;
			}

			std::cout << std::endl;
			std::cout << std::endl;
			std::cout << std::endl;
			dumpMatrix(z[0], cols, rows);
			dumpMatrix(z[0], SIZE, rows);
			calculate2(z[1], z[0], cols, rows);
			Comm_Info info;
			setSyncInfo(info, z[1], cols, rows);
			communicateOverlap(info);
			dumpMatrix(z[1], cols, rows);
			dumpMatrix(z[1], SIZE, rows);
			reconstructMatrix(z[1], SIZE, SIZE);
			dumpMatrix(z[1], cols, rows);
			dumpMatrix(z[1], SIZE, rows);

			for (auto i = 1; i <= SIZE; i++)
			{
				for (auto j = 0; j < SIZE; j++)
				{
					if (z[1][i][j] < 10)
					{
						std::cout << z[1][i][j] << "    ";
					}
					else if (z[1][i][j] < 100)
					{
						std::cout << z[1][i][j] << "   ";
					}
					else if (z[1][i][j] < 1000)
					{
						std::cout << z[1][i][j] << "  ";
					}
					else
					{
						std::cout << z[1][i][j] << " ";
					}
				}
				std::cout << std::endl;
			}

		} // end if
		else
		{
			dumpMatrix(z[0], cols, rows);
			calculate2(z[1], z[0], cols, rows);
			Comm_Info info;
			setSyncInfo(info, z[1], cols, rows);
			communicateOverlap(info);
			dumpMatrix(z[1], cols, rows);
			reconstructMatrix(z[1], cols, rows);
			dumpMatrix(z[1], cols, rows);
		} // end else

		MPI_Finalize();

		deallocate(z, SIZE);
	}


	void calculate2(double** t_0, double** t_1, const int COLS, const int ROWS)
	{
		for (auto i = 1; i <= COLS; i++)
		{
			for (auto j = 0; j < ROWS; j++)
			{
				// manager has left boundary
				if (!my_rank)
				{
					if (i == 1 || !j || j + 1 == ROWS)
					{
						t_0[i][j] = 1.0;
					} // end if
					else
					{
						t_0[i][j] = (t_1[i + 1][j] + t_1[i - 1][j] + t_1[i][j + 1] + t_1[i][j - 1]);
					} // end else
				} // end if
				else if (my_rank + 1 == mpi_size) // right-most worker has right boundary
				{
					if (!j || i == COLS || j + 1 == ROWS)
					{
						t_0[i][j] = 2.0;
					} // end if
					else
					{
						t_0[i][j] = (t_1[i + 1][j] + t_1[i - 1][j] + t_1[i][j + 1] + t_1[i][j - 1]);
					} // end else
				} // end elif
				else // all other ranks have no right and left boundaries
				{
					if (!j || j + 1 == ROWS)
					{
						t_0[i][j] = 3.0;
					} // end if
					else
					{
						t_0[i][j] = (t_1[i + 1][j] + t_1[i - 1][j] + t_1[i][j + 1] + t_1[i][j - 1]);
					} // end else
				} // end else
			} // end for j
		} // end for i
	} // end method calculate(3)


	void dumpMatrix(double** t_0, const int COLS, const int ROWS)
	{
		std::stringstream ss;

		ss << "mat_dump_" << my_rank << ".txt";

		std::string fname = ss.str();
		std::ofstream mat_dump(fname, std::ios::out | std::ios::app);

		for (auto i = 0; i < COLS + 2; i++)
		{
			for (auto j = 0; j < ROWS; j++)
			{
				if (t_0[i][j] < 10)
				{
					mat_dump << t_0[i][j] << "    ";
				}
				else if (t_0[i][j] < 100)
				{
					mat_dump << t_0[i][j] << "   ";
				}
				else if (t_0[i][j] < 1000)
				{
					mat_dump << t_0[i][j] << "  ";
				}
				else
				{
					mat_dump << t_0[i][j] << " ";
				}
			}
			mat_dump << "\n";
		}

		mat_dump << "\n";
		mat_dump << "\n";
		mat_dump << "\n";
		mat_dump.close();
	}

#endif

#pragma endregion


