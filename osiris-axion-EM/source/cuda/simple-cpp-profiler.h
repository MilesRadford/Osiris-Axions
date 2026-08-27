    #ifndef GHETTO_PROFILER_H
    #define GHETTO_PROFILER_H
	
#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
	#define WINDOWS
#endif

#ifdef WINDOWS
	//#include "stdafx.h"
	#include <time.h>
	#include <windows.h>
#else
	#include <sys/time.h>
#endif
	
  #include <iosfwd>
	#include <ostream>
  #include <iostream>
  #include <fstream>
	#include <iomanip>
	#include <stack>
  #include <unordered_map>
  #include <vector>
	#include <sstream>


#ifdef WINDOWS
	#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
		#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
	#else
		#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
	#endif

	struct timezone
	{
	  int  tz_minuteswest; /* minutes W of Greenwich */
	  int  tz_dsttime;     /* type of dst correction */
	};

#endif

//using namespace std;
	/*
		TODO: This is code is close but not quite THREAD SAFE. But since MPI really doesnt do multi-threaded executuion
				we will be ok.. but I will stick a barrier block on a few of these methods sometime and then it will be ok.
				
		TODO: Maybe I if get bored, replace string hashes with enumerated types... but I like strings as hashtable keys.. 
				makes me feel like it's Java.
		
		This works better the more tree like the calls are... in partickualar, each function is 'frozen' at the call stack depth 
		it was first seen.. so if 2 function at level 2 is later called by a funciton at level 10 then the time it takes will
		be charged as if it was occuring at level 2.. 
		
		So this code really deigned for the first few (3-4) call stack levels from main() for high level algorithm
		decisions.
	*/
	

	/* prolly a micro-optimzation.. but it's easy.. add an integer event id (in addition to the string one). The string one is nice because you
		can just put it anywhere you want to measure and that's it. But is had a bunch of string manipulations in it that can start to effect speed
		of the code being tested if it's a tught loop... so for cases where you would leik fast-as-possible or for cases where the timting will be left permanently in the 
		code, add another method using ineger enumberated types.. this is slightly less convienent then the string calss as you have to register 
		any new events here..
	
	enum ghetto_timer_events {
		CUDA_ee_collsions,
		CPU_ee_collsions,
		ee_collsion_step,
		ee_collsions_RK4_calculations,
		ee_collsion_per_cell,
		Actions_Adv,
		Implicit_Actions_Adv1,
		Implicit_Actions_Adv2,
		Implicit_Actions_ImpEx,
		Implicit_Actions_ImpEy,
		Implicit_Actions_ImpEz,
		on_da_CLEAN_down_LOW_PC,
		E_Field_Meth_Implicit_E_Field,
		Node_ImplicitE_Communication,
		PE_Neighbor_Communications,
		Rk2_Implicit_p1,
		Rk2_Implicit_p2,
		ElectricField_Exyz1D,
		ElectricField_Exyz,
		MagneticField_Bxyz,
		Spatial_Advection_Axy,
		MaxwellEq_Faraday,
		Main_Loop,
		ee_collsion_test_overall
	};

	const char* ghett_event_labels[] = {
		"CUDA e*e collsions",
		"CPU e*e collsions",
		"e*e collsion step",
		"e*e collsions RK4 calculations",
		"e*e collsion per cell",
		"Actions::Adv",
		"Implicit_Actions::Adv1",
		"Implicit_Actions::Adv2",
		"Implicit_Actions::ImpEx",
		"Implicit_Actions::ImpEy",
		"Implicit_Actions::ImpEz",
		"on_da_CLEAN_down_LOW_PC",
		"E_Field_Meth::Implicit_E_Field",
		"Node_ImplicitE_Communications",
		"PE::Neighbor_Communications"
		"Rk2 Implicit p1",
		"Rk2 Implicit p2",
		"ElectricField::Exyz1D",
		"ElectricField::Exyz",
		"MagneticField::Bxyz",
		"Spatial_Advection::Axy",
		"MaxwellEq::Faraday",
		"Main Loop",
		"e*e collsion test overall"
	}; */
	
	// declaration of simple Log class (writing to disk)
	class Log	{
		private:
			int my_rank;
			
		public:
			std::ofstream out;
			//std::ostream logger;
			Log();
			Log(int my_rank);
			~Log();
			void init(int rank);
			void init(const char* name, int rank);
			void l(std::string message);
			template <class T>
				void dump_buffer(T* buffer, int len, std::string prefix, int rank);
	};
 	
 	// forward declaration of static LOG instance
	std::ostream& LOG();
  

  // Aggregate information kept for each distinct timer
	class timer_residue {
		public: 
    double dt;
		double accumlated_time;
		unsigned long times_used;
		double maxx;
		double minn;
		size_t call_stack_depth;
		timer_residue* parent;
		int parent_ordinal;
		std::vector<timer_residue*> children();
		std::string name;
		int ordinal;
		timer_residue(): dt(0.0), accumlated_time(0.0), times_used(0), maxx(0.0), minn(0.0), call_stack_depth(0),ordinal(0) {}
	};
	
	// this is a quick forward declaration......
	class ScopeTimer;	
	
	
	/*
		Singleton class that control and colates all the timer results.
		There is - and only can be - One TimerOracle.
	*/
	class TimerOracle {
		public:
			typedef std::unordered_map<std::string, timer_residue*> TimerData;
			TimerData __timer_data;
			std::stack<timer_residue*> call_state_stack;
			int number_nodes_added;
			
			bool enabled;
			bool write_output;
			
			std::string current_context;
			
			// makes our Oracle a singularity.. I mean a Singleton.
			static TimerOracle* summon();

			timer_residue* hearld_birth_stack (std::string name) {
				if(call_state_stack.empty()) {
					std::string composite_key;
					timer_residue* timer_entry;
					if(__timer_data.count(composite_key)) {
						// it was in hashtbale, so retrive it..
						timer_entry = __timer_data[composite_key];
					} else {
						// not in hashtable.. make new object and put in hashtable.
						timer_entry = new timer_residue();
						__timer_data[composite_key] = timer_entry;
						timer_entry->name = (name);
						timer_entry->ordinal = number_nodes_added++;
					}
					timer_entry->call_stack_depth = 0;
					timer_entry->parent = NULL;
					timer_entry->parent_ordinal = -1;
					timer_entry->name = (name);
					call_state_stack.push(timer_entry);
					return timer_entry;
				} else {
					timer_residue* &current_parent = call_state_stack.top();
					// does he know that I am his child?
					// make a composite key (slowish.. fix later)
					std::string composite_key = current_parent->name + name;

					timer_residue* timer_entry;
					if(__timer_data.count(composite_key)) {
						// it was in hashtbale, so retrive it..
						timer_entry = __timer_data[composite_key];
					} else {
						// not in hashtable. make new object and put in hashtable.
						timer_entry = new timer_residue();
						__timer_data[composite_key] = timer_entry;
						timer_entry->name = (name);
						timer_entry->ordinal = number_nodes_added++;
					}
					
					timer_entry->call_stack_depth = call_state_stack.size();
					timer_entry->parent = current_parent;
					timer_entry->parent_ordinal = current_parent->ordinal;
					call_state_stack.push(timer_entry);
					return timer_entry;
				}
			}
		
		void sacrfice_offering_stack(double dt, timer_residue* timer_entry) {
			//if(enabled == false) return;
			if(timer_entry == NULL) return;
			
			(*timer_entry).dt = dt;
			(*timer_entry).accumlated_time += dt;
			(*timer_entry).times_used++;
      
      call_state_stack.pop();
		}
		

		timer_residue* hearld_birth(std::string name) {

			// get the "timer_residue" out of the hashtable. If it is not found,
			// 		make a new entry in the hashtable.
			timer_residue* timer_entry;
			if(__timer_data.count(name)) {
				// it was in hashtable, so retrive it..
				timer_entry = __timer_data[name];
			} else {
				// not in hastbale.. make new object and put in hashtable.
				timer_entry = new timer_residue();
				__timer_data[name] = timer_entry;
				timer_entry->name = (name);
			}
			
			// this is start of a stack implementation.. for now, at least get relative stack depth.
			// 		we can make it as fancy was we want later if there's a need.
			return timer_entry;
		}
		
		// fetcher functions for the stored timer residues...
		timer_residue* get_result(const char * name) {
			return get_result(std::string(name));
		}


		timer_residue* get_result(std::string name) {
			std::string composite_key;
			if(call_state_stack.empty()) {
				composite_key.assign(name);
			} else {
				timer_residue* &current_parent = call_state_stack.top();
				composite_key.assign(current_parent->name + name);
			}

			if(__timer_data.count(composite_key)) {
				return __timer_data[composite_key];
			}

			return NULL;
		}
		
		inline bool is_enabled() {
			return this->enabled;
		}
		
		void global_disable_timing() {
			this->enabled = false;
		}
		void global_enable_timing() {
			this->enabled = true;
		}
		
		void dump_nvpair_new_format(double time, int time_step) {
			if(write_output)
				LOG() << "<<profile_timings={";
			TimerData::const_iterator current = __timer_data.begin(); 
			TimerData::const_iterator end = __timer_data.end(); int i=0;
			while(current != end) {
				if(write_output) {
					if(i > 0) {
						LOG() << ", '" << ((current)->second)->name << "':(" << ((current)->second)->accumlated_time << "," << ((current)->second)->times_used << "," << ((current)->second)->call_stack_depth << ","  << ((current)->second)->ordinal <<  ","  << ((current)->second)->parent_ordinal << ")";
					} else {
						LOG() << "'" << ((current)->second)->name << "':(" << ((current)->second)->accumlated_time << "," << ((current)->second)->times_used << "," << ((current)->second)->call_stack_depth << ","  << ((current)->second)->ordinal <<  ","  << ((current)->second)->parent_ordinal << ")";
					}
				}
				
				// zero this out after flushing just keep precision on the floating point numbers.
				((current)->second)->accumlated_time = 0.0;
				((current)->second)->times_used = 0;
				current++; i++;
			}
			
			if(write_output) {
				LOG() << "}" << std::endl;
				LOG() << "<<simulation_times={'time':" << time << ", 'time_step':" <<  time_step << "}" << std::endl;
			}
			
		}
		
		void change_context( std::string new_context, double time, int time_step, int num_particles_this_node, double num_particles__sum_all_nodes__all__species, int num_nodes, std::string old_params, std::string old_context_params) {
			
			this->dump_nvpair(time, time_step, num_particles_this_node, num_particles__sum_all_nodes__all__species, num_nodes, old_params, old_context_params);
			this->current_context = new_context;
		}
		
		void dump_nvpair(double time, int time_step, int num_particles_this_node, double num_particles__sum_all_nodes__all__species, int num_nodes, std::string params, std::string context_params) {
						
			if(write_output) {
				LOG() << "???Dataset=['context','context_params','Global_Timings','time','num_par', 'params']"<< std::endl;
				LOG() << "{ 'context': '" << this->current_context << "' }";
				LOG() << std::endl << context_params;
				LOG() << std::endl << "{";
			}
			TimerData::const_iterator current = __timer_data.begin(); 
			TimerData::const_iterator end = __timer_data.end(); int i=0;
			while(current != end) {
				if(write_output) {
					if(i > 0) {
						LOG() << ", '" << ((current)->second)->name << "':(" << ((current)->second)->accumlated_time << "," << ((current)->second)->times_used << "," << ((current)->second)->call_stack_depth << ","  << ((current)->second)->ordinal <<  ","  << ((current)->second)->parent_ordinal << ")";
					} else {
						LOG() << "'" << ((current)->second)->name << "':(" << ((current)->second)->accumlated_time << "," << ((current)->second)->times_used << "," << ((current)->second)->call_stack_depth << ","  << ((current)->second)->ordinal <<  ","  << ((current)->second)->parent_ordinal << ")";
					}
				}
				
				// zero this out after flushing just keep precision on the floating point numbers.
				((current)->second)->accumlated_time = 0.0;
				((current)->second)->times_used = 0;
				current++; i++;
			}
			if(write_output) {
				LOG() << "}"; 
				LOG() << std::endl << "{'time':" << time << ", 'time_step':" <<  time_step << "}";
				LOG() << std::endl << "{'num_par':" << num_particles_this_node << ", 'total':" <<  num_particles__sum_all_nodes__all__species << ",'num_nodes':" << num_nodes << "}";
				LOG() << std::endl << params << std::endl;
				LOG() << "???ItIsFinished"<<std::endl;
			}
			
		}
    
    void force_parent(timer_residue* child, timer_residue* parent) {
      std::string composite_key = parent->name + child->name;
      
      child->parent = parent;
      child->parent_ordinal = parent->ordinal;
      
      __timer_data[composite_key] = child;
      __timer_data.erase(child->name);      
    }
	
		/*		for 'jagged (non-stack structured) timers which we dont need yet (or maybe ever)
		TimerAgent* register_timer(string name) {
			__registered_timers[name] = new TimerAgent();
			return __registered_timers[name];
		}
		TimerAgent* register_timer(string context, string name) {
			__registered_timers[name] = new TimerAgent();
			return __registered_timers[name];
		}
		TimerAgent* get_timer(string name) {
			return  __registered_timers[name];
		}*/
			
			
		TimerOracle() : current_context("Default") {number_nodes_added=0; enabled=false; write_output = false; };
		TimerOracle(TimerOracle const&);				// prevent copying of the Oracle.
		TimerOracle& operator=(TimerOracle const&);		// prevent addigning the Oracle to different avatars.
		static TimerOracle* m_pInstance;
		~TimerOracle();
	};
	
	
	/*
		The little stop watches we plant around the code to keep track of performance. Note, this timer is meant
		for timing Scopes... it starts as soon as this object is created... and stops when the object goes out of scope and
		destructor is called.
	*/
	class ScopeTimer {
		private:
			timeval t_start;
			timeval t_end;
			timer_residue* pointer_to_where_my_soul_will_go_when_I_die;
		public:
			ScopeTimer(const char * name) {
					if( TimerOracle::summon()->enabled)
						pointer_to_where_my_soul_will_go_when_I_die = TimerOracle::summon()->hearld_birth_stack(name);
					else
						pointer_to_where_my_soul_will_go_when_I_die = NULL;
					gettimeofday(&t_start, NULL); 
			}
			ScopeTimer(std::string name) { 
				if( TimerOracle::summon()->enabled)
					pointer_to_where_my_soul_will_go_when_I_die = TimerOracle::summon()->hearld_birth_stack(name);
					else
						pointer_to_where_my_soul_will_go_when_I_die = NULL;
				
				gettimeofday(&t_start, NULL); 
			}
		~ScopeTimer() {
				std::cout << "ScopeTimer Dieid" << std::endl;
				if( pointer_to_where_my_soul_will_go_when_I_die != NULL) {
						std::cout << "ScopeTimer Dieid with hinirs" << std::endl;
						gettimeofday(&t_end, NULL);
            // get time and convert it to milliseconds (tv_sec gives seconds, tv_usec give microseconds)
						double dt = (t_end.tv_sec - t_start.tv_sec) * 1000.0 + (t_end.tv_usec - t_start.tv_usec) / 1000.0; 
						TimerOracle::summon()->sacrfice_offering_stack(dt, pointer_to_where_my_soul_will_go_when_I_die); 
				}
		}
		
		
		#ifdef WINDOWS
			int gettimeofday(struct timeval *tv, struct timezone *tz)
			{
			// Define a structure to receive the current Windows filetime
			  FILETIME ft;
			 
			// Initialize the present time to 0 and the timezone to UTC
			  unsigned __int64 tmpres = 0;
			  static int tzflag = 0;
			 
			  if (NULL != tv) {
			    GetSystemTimeAsFileTime(&ft);
			 
				// The GetSystemTimeAsFileTime returns the number of 100 nanosecond 
				// intervals since Jan 1, 1601 in a structure. Copy the high bits to 
				// the 64 bit tmpres, shift it left by 32 then or in the low 32 bits.
			    tmpres |= ft.dwHighDateTime;
			    tmpres <<= 32;
			    tmpres |= ft.dwLowDateTime;
			 
				// Convert to microseconds by dividing by 10
			    tmpres /= 10;
			 
				// The Unix epoch starts on Jan 1 1970.  Need to subtract the difference 
				// in seconds from Jan 1 1601.
			    tmpres -= DELTA_EPOCH_IN_MICROSECS;
			 
				// Finally change microseconds to seconds and place in the seconds value. 
				// The modulus picks up the microseconds.
			    tv->tv_sec = (long)(tmpres / 1000000UL);
			    tv->tv_usec = (long)(tmpres % 1000000UL);
			  }
			 
			  if (NULL != tz) {
			    if(!tzflag) {
					_tzset();
					tzflag++;
				}
				// Adjust for the timezone west of Greenwich
				//tz->tz_minuteswest = _timezone / 60;
			    //tz->tz_dsttime = _daylight;
			  }
			  return 0;
			}
		#endif
	};
	

	class StopWatchTimer {
		private:
			timeval t_start;
			timeval t_end;
			timer_residue* pointer_to_where_my_soul_will_go_when_I_die;
			bool clock_paused;
			
		public:
			StopWatchTimer(std::string name) : clock_paused(false) { 
				pointer_to_where_my_soul_will_go_when_I_die = TimerOracle::summon()->hearld_birth_stack(name);
				gettimeofday(&t_start, NULL); 
			}
			StopWatchTimer(const char * name) : clock_paused(false) { 
				pointer_to_where_my_soul_will_go_when_I_die = TimerOracle::summon()->hearld_birth_stack(name);
				gettimeofday(&t_start, NULL); 
			}

						
			void force_parent(timer_residue* parent_pointer) {
				TimerOracle::summon()->force_parent(pointer_to_where_my_soul_will_go_when_I_die,parent_pointer); 
			  
			}
			
			void pause_clock() {
			  // get time and convert it to milliseconds (tv_sec gives seconds, tv_usec give microseconds)
			  gettimeofday(&t_end, NULL);
			  double dt = (t_end.tv_sec - t_start.tv_sec) * 1000.0 + (t_end.tv_usec - t_start.tv_usec) / 1000.0;
			  if( !TimerOracle::summon()->enabled) return;
			  (*pointer_to_where_my_soul_will_go_when_I_die).dt = dt;
			  (*pointer_to_where_my_soul_will_go_when_I_die).accumlated_time += dt;
			  
			  // Note: we specifically do not updatte the discrete interval count in this case.
			  //(*pointer_to_where_my_soul_will_go_when_I_die).times_used++;
			  clock_paused = true;
			}
			void unpause_clock() {
				clock_paused = false;
				gettimeofday(&t_start, NULL);
			}
			void zero_clock() {
			  (*pointer_to_where_my_soul_will_go_when_I_die).dt = 0;
			  (*pointer_to_where_my_soul_will_go_when_I_die).accumlated_time = 0;
			  gettimeofday(&t_start, NULL); 
				std::cout << " zero clock sucks.." << std::endl;
			  exit(1);
			}
			
			void mark_interval() {
			  
			  // get time and convert it to milliseconds (tv_sec gives seconds, tv_usec give microseconds)
			  gettimeofday(&t_end, NULL);
			  double dt = (t_end.tv_sec - t_start.tv_sec) * 1000.0 + (t_end.tv_usec - t_start.tv_usec) / 1000.0;
			  if( !TimerOracle::summon()->enabled) return;
			  
			  if(!clock_paused) {
				(*pointer_to_where_my_soul_will_go_when_I_die).dt = dt;
				(*pointer_to_where_my_soul_will_go_when_I_die).accumlated_time += dt;
			  }
			  
			  (*pointer_to_where_my_soul_will_go_when_I_die).times_used++;
			  
			  t_start.tv_sec = t_end.tv_sec;
			  t_start.tv_usec = t_end.tv_usec;
			  //LOG() << "IN MARK TIM AT LOOP END! ITS MY TRIP NUMBER " << (*pointer_to_where_my_soul_will_go_when_I_die).times_used << " tot time: " << (*pointer_to_where_my_soul_will_go_when_I_die).accumlated_time << endl << endl;
			}
			
			void mark_interval(int number_of_intervals) {
			  
			  // get time and convert it to milliseconds (tv_sec gives seconds, tv_usec give microseconds)
			  gettimeofday(&t_end, NULL);
			  double dt = (t_end.tv_sec - t_start.tv_sec) * 1000.0 + (t_end.tv_usec - t_start.tv_usec) / 1000.0;
			  if( !TimerOracle::summon()->enabled) return;
			  
			  if(!clock_paused) {
				(*pointer_to_where_my_soul_will_go_when_I_die).dt = dt;
				(*pointer_to_where_my_soul_will_go_when_I_die).accumlated_time += dt;
			  }
			  
			  (*pointer_to_where_my_soul_will_go_when_I_die).times_used += number_of_intervals;
			  t_start.tv_sec = t_end.tv_sec;
			  t_start.tv_usec = t_end.tv_usec;
			  //LOG() << "IN MARK TIM AT LOOP END! ITS MY TRIP NUMBER " << (*pointer_to_where_my_soul_will_go_when_I_die).times_used << " tot time: " << (*pointer_to_where_my_soul_will_go_when_I_die).accumlated_time << endl << endl;
			}
			
			int get_current_number_intervals_measured() {
				return (*pointer_to_where_my_soul_will_go_when_I_die).times_used;
			}
			
			~StopWatchTimer() {
				gettimeofday(&t_end, NULL);
				// get time and convert it to milliseconds (tv_sec gives seconds, tv_usec give microseconds)
				double dt = (t_end.tv_sec - t_start.tv_sec) * 1000.0 + (t_end.tv_usec - t_start.tv_usec) / 1000.0; 
				TimerOracle::summon()->sacrfice_offering_stack(dt, pointer_to_where_my_soul_will_go_when_I_die); 
			}
		
		
		#ifdef WINDOWS
			int gettimeofday(struct timeval *tv, struct timezone *tz) {
			// Define a structure to receive the current Windows filetime
			  FILETIME ft;
			 
			// Initialize the present time to 0 and the timezone to UTC
			  unsigned __int64 tmpres = 0;
			  static int tzflag = 0;
			 
			  if (NULL != tv) {
				GetSystemTimeAsFileTime(&ft);
			 
				// The GetSystemTimeAsFileTime returns the number of 100 nanosecond 
				// intervals since Jan 1, 1601 in a structure. Copy the high bits to 
				// the 64 bit tmpres, shift it left by 32 then or in the low 32 bits.
				tmpres |= ft.dwHighDateTime;
				tmpres <<= 32;
				tmpres |= ft.dwLowDateTime;
			 
				// Convert to microseconds by dividing by 10
				tmpres /= 10;
			 
				// The Unix epoch starts on Jan 1 1970.  Need to subtract the difference 
				// in seconds from Jan 1 1601.
				tmpres -= DELTA_EPOCH_IN_MICROSECS;
			 
				// Finally change microseconds to seconds and place in the seconds value. 
				// The modulus picks up the microseconds.
				tv->tv_sec = (long)(tmpres / 1000000UL);
				tv->tv_usec = (long)(tmpres % 1000000UL);
			  }
			 
			  if (NULL != tz) {
				if(!tzflag) {
					_tzset();
					tzflag++;
				}
				// Adjust for the timezone west of Greenwich
				//tz->tz_minuteswest = _timezone / 60;
				//tz->tz_dsttime = _daylight;
			  }
			  return 0;
			}
		#endif
  };
	
	// forward declarations
	Log& LOG_OBJECT();
	std::ostream& LOG();

	#endif