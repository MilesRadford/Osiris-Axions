//using namespace std;

#include "simple-cpp-profiler.h"


// Simple Singleton pattern (TODO: make a real singleton for fun)
TimerOracle* TimerOracle::m_pInstance = NULL; 

TimerOracle* TimerOracle::summon() {
	if (!m_pInstance)   // Only allow one instance of class to be generated.
		m_pInstance = new TimerOracle();
	return m_pInstance;
}

// The following is the defintio
//struct log_nullstream {};nl of a null stream for disabling the logging..
//template <typename T>
//log_nullstream& operator << (log_nullstream& s, T const & ) { return s; }
//log_nullstream& operator << (log_nullstream& s, std::ostream &(std::ostream& )) { return s; }
//
//static log_nullstream nulllogstream;
	
	
	Log::Log() {
		my_rank = -1;
	}
	
	Log::~Log() {
		//cout << "WHo killed my LOG?!?!?!" << endl;
		//out.close();
	}
	 
	Log::Log(int rank) {
		//if(disabled) {
		//	logger = nulllogstream;
		//} else {
			my_rank = rank;
			std::stringstream fname;
			fname << "log_" << rank << ".txt";
			out.open(fname.str().c_str());
			//logger = out;
		//}
	}
	void Log::init(int rank) {
		//if(disabled) {
		//	logger = nulllogstream;
		//} else {
			my_rank = rank;
			std::stringstream fname;
			fname << "log_" << rank << ".txt";
			out.open(fname.str().c_str());
			//logger = out;
		//}
	}
	
	void Log::init(const char* name, int rank) {
		//if(disabled) { 
		// logger = nulllogstream;
		//} else {
			my_rank = rank;
			std::stringstream fname;
			fname << "log_" << name << "_" << rank << ".txt";
			out.open(fname.str().c_str());
			//logger = out;
		//}
	}
	void Log::l(std::string s) {
		out << s << std::endl;
		out.flush();
	}
	
    Log& LOG_OBJECT() {
        static Log log;
        return log;
    }
	std::ostream& LOG() {
		return LOG_OBJECT().out;
	}
	
    Log& LOG1_OBJECT() {
        static Log log;
        return log;
    }
	std::ostream& LOG1() {
		return LOG1_OBJECT().out;
	}
	
    Log& LOG2_OBJECT() {
        static Log log;
        return log;
    }
	std::ostream& LOG2() {
		return LOG2_OBJECT().out;
	}