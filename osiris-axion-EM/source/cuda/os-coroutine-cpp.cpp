/*
    TODO: 
        1. check to make sure no memory leaks are occuring.
        2. SIM1_ITER, SIM2_ITER should be function pointers..
            I am jsut beign lazy for now about getting Fortran to pass the correct pointers.
*/

#ifdef ENABLE_COROUTINE

#include <iostream>
#include "boost/coroutine2/all.hpp"
#include <string>


using namespace std;
typedef boost::coroutines2::coroutine<string>   coro_t;
typedef boost::coroutines2::detail::pull_coroutine<string>::iterator   coro_t_itt;

extern "C" void SIM1_ITER();
extern "C" void SIM2_ITER();

struct SimulationContextManager {
  coro_t::pull_type* source__1;
  coro_t::pull_type* source__2;
  coro_t::push_type* sink__1;
  coro_t::push_type* sink__2;

  coro_t_itt it_1;
  coro_t_itt it_2;
  int current_context = 0;
  string last_returned_context;
};

// make this a global for now.. easy enough to make relevent local variables if ever required.
SimulationContextManager __contextManager;

extern "C" void init_cosimulation() {
  cout << "----------------" << endl;
  cout << "HI.. This is init_cosimulation" << endl;
  cout << "----------------" << endl;
}

extern "C" void start_sim_itteration() {
    // coro_t::pull_type source(
    //   [&](coro_t::push_type& sink){
    //       SIM1_ITER();
    //       //string first = string("A");
    //       //sink(first);
    //       //for(int i=0;i<8;++i){
    //       //    first = first + string("A");
    //       //    sink(first);
    //       //}
    //   });
}

extern "C" void start_sim_itteration_cpp() {

  // create the coroutine for SIM1
  __contextManager.source__1 = new coro_t::pull_type( [&](coro_t::push_type& sink){
      string temp;

      // save a pointer to the 'sink__1' function. This is the function that is called
      //    to pause the execution and return to the calling context
      //    (in our case, SIM1 will call this function to pause itself)
      __contextManager.sink__1 = &sink;

      // The following is a throwaway 'primer' to immediatly pause+return to the calling context 
      //    as this coroutine is created so we can synchronize.
      sink(string("__JUNK"));
      SIM1_ITER();
    });

  // create the coroutine for SIM2
  __contextManager.source__2 = new coro_t::pull_type( [&](coro_t::push_type& sink){
      string temp;

      // save a pointer to the 'sink__2' function. This is the function that is called
      //    to pause the execution and return to the calling context
      //    (in our case, SIM2 will call this function to pause itself)
      __contextManager.sink__2 = &sink;

      // The following is a throwaway 'primer' to immediatly pause+return to the calling context 
      //    as this coroutine is created so we can synchronize.
      sink(string("__JUNK"));
      SIM2_ITER();
    });

  // create C++ iterators for advancing though the coroutine contexts
  __contextManager.it_1 = begin( *(__contextManager.source__1) );
  __contextManager.it_2 = begin( *(__contextManager.source__2) );

  // initalize our internal state variable identifing the most recent coroutine context run.
  __contextManager.current_context = -1;

}

void do_iter(coro_t_itt *it, coro_t::pull_type* source ) {

  // If the coroutine simply finished (rather then requesting a pause),
  //  then report a context id of '__END'
  if( (*it) == end(*source) ) {
    __contextManager.last_returned_context = "__END";
    return;
  }

  // Run this coroutine starting immedialy after the pause point.
  (*it)++;

  // If the coroutine simply finished (rather then requesting a pause),
  //  then report a context id of '__END'
  if( (*it) == end(*source) ) {
    __contextManager.last_returned_context = "__END";
  } else {
    // The code has paused itself! Save the context id the code paused itself in.
    __contextManager.last_returned_context = *(*it);
  }

}

extern "C" void execute_sim1_iter_cpp( ) {
  __contextManager.current_context = 1;
  do_iter(&(__contextManager.it_1), __contextManager.source__1);
}

extern "C" void execute_sim2_iter_cpp( ) {
  __contextManager.current_context = 2;
  do_iter(&(__contextManager.it_2), __contextManager.source__2);
}


// This is the primary routine called within the Fortran code that
//  causes execution to paused.
extern "C" void handoff_context_cpp(const char *path) {
   
  if(__contextManager.current_context == 1) {
    (*(__contextManager.sink__1))( string(path) );
  }
  if(__contextManager.current_context == 2) {
    (*(__contextManager.sink__2))( string(path) );
  }  
}


/*
    Copies the current context id ('last_returned_context') to a Fortran string.
      Fortran strings are fixed length and are ppadded out to this length 
      with ' ' (that is, space) characters.
*/
extern "C" void test_str_cpp(char *path) {
  std::strcpy (path,  __contextManager.last_returned_context.c_str());
  for(int i=__contextManager.last_returned_context.length();i<50;i++) {
    path[i] = ' ';
  }
}

#endif


























#if 0
// This is test BOOST code that went into making this

#include <iostream>
//#include <boost\coroutine2\all.hpp>
#include <boost/coroutine2/all.hpp>
#include <unordered_map>
#include <map>
#include <string>

using namespace std;

void test1() {

  /*
  typedef boost::coroutines2::coroutine<int>   coro_t;
  typedef boost::coroutines2::detail::pull_coroutine<int>::iterator   coro_t_itt;

  coro_t::pull_type source(
      [&](coro_t::push_type& sink){
          int first=1,second=1;
          sink(first);
          sink(second);
          for(int i=0;i<8;++i){
              int third=first+second;
              first=second;
              second=third;
              sink(third);
          }
      }); 
  */

  typedef boost::coroutines2::coroutine<string>   coro_t;
  typedef boost::coroutines2::detail::pull_coroutine<string>::iterator   coro_t_itt;

  coro_t::pull_type source(
      [&](coro_t::push_type& sink){
          string first = string("A");
          sink(first);
          for(int i=0;i<8;++i){
              first = first + string("A");
              sink(first);
          }
      });

  std::map<string, coro_t_itt> contexts;

  contexts[string("CONTEXT1")] = begin(source);
  //for(auto i:source)
  //  std::cout << i <<  " ";
  //coro_t_itt it = begin(source);
  //auto nx = std::next(it);

  coro_t_itt it = contexts[string("CONTEXT1")];

  cout << *(contexts[string("CONTEXT1")]) << "\n";
  (contexts[string("CONTEXT1")])++;

  cout << *(contexts[string("CONTEXT1")]) << "\n";
  (contexts[string("CONTEXT1")])++;

  cout << *(contexts[string("CONTEXT1")]) << "\n";
  (contexts[string("CONTEXT1")])++;

  cout << *(contexts[string("CONTEXT1")]) << "\n";
  (contexts[string("CONTEXT1")])++;

  cout << *it << "\n";
  it++;
  cout << *it << "\n";

  it++;
  cout << *it << "\n";

  it++;
  cout << *it << "\n";

  it++;
  cout << *it << "\n";

  (contexts[string("CONTEXT1")])++;
  if((contexts[string("CONTEXT1")]) == end(source)) cout << "end 2!" << "\n";

  cout << "sdfsdf";
}

void test2() {
  typedef boost::coroutines2::coroutine<string>   coro_t;
  typedef boost::coroutines2::detail::pull_coroutine<string>::iterator   coro_t_itt;

  coro_t::pull_type* source__1;

  coro_t_itt *it_1;

  source__1 = new coro_t::pull_type([&](coro_t::push_type& sink){
          string first = string("A");
          sink(first);
          for(int i=0;i<8;++i){
              first = first + string("A");
              sink(first);
          }
      });
  (it_1) = &(begin( *source__1 ));
  (*it_1)++;
  cout << *(*it_1) << endl;
  (*it_1)++;
  cout << *(*it_1) << endl;
  (*it_1)++;
  cout << *(*it_1) << endl;
  // auto it_11 = begin( *source__1 );
  // it_11++;
  // cout << *it_11 << endl;
  // it_11++;
  // cout << *it_11 << endl;
  // it_11++;
  // cout << *it_11 << endl;
  // cout << *it_11 << endl;
  // cout << *it_11 << endl;
}

void main() {
  //test1();
  test2();

}

#endif

#if 0
// This is test CMAKE file (CmakeLists.txt)
cmake_minimum_required(VERSION 2.8)

project (boostme)


include_directories( "C:/lab_code/boost_1_70_0_b1_rc3")

add_executable(boostme "test.cpp")
TARGET_LINK_LIBRARIES(boostme "C:/lab_code/boost_1_70_0_b1_rc3/lib64-msvc-14.1/libboost_context-vc141-mt-gd-x64-1_70.lib")

#endif