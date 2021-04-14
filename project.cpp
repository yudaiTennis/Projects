#include <stdio.h> 
#include <stdlib.h>  
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include <math.h>
#include <fstream>
#include <queue>
#include <utility>
#include <algorithm>
#include <numeric>

using namespace std;

#define process_id "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

class Process;
void print_queue(vector<Process> queue);
void print_queue(priority_queue<Process> queue);
vector<double> fcfs(vector<Process> processes, int n, int context_time);
vector<double> sjf(vector<Process> processes, int n, int context_time, int lambda);
vector<double> srt(vector<Process> processes, int n, int context_time, int lambda);
vector<double> rr( vector<Process> processes, int n, int timeslice, bool rr_add, int context_time); 

class Process
{
public:
    char id;
    int arrival_time;
    bool firstTime;
    int init_number_of_bursts;
    vector<int> init_cpu;
    int working_burst;
    int wait_time;
    vector<int> tau;
    vector<int>  cpu_bursts;
    vector<int> io_bursts;
    bool preempt;

    Process(int i, int a_time, int n_bursts, vector<int> cpu, vector<int> io, double lambda, double alpha)
    {
        id = process_id[i];
        arrival_time = a_time;
        init_number_of_bursts = n_bursts;
        working_burst = 0;
        init_cpu = cpu;
        cpu_bursts = cpu;
        io_bursts = io;
        firstTime = true;
        wait_time = 0;
        set_tau(lambda, alpha);
        preempt = false;
    }
    Process(void)
    {
        id = '.';
        arrival_time = 0;
        init_number_of_bursts = 0;
        working_burst = 0;
        wait_time = 0;
    }

    void startIO(int t, vector<Process> queue) {
        arrival_time += io_bursts[working_burst];
        if (t < 1000) {cout << "time " << t << "ms: Process "<< id << " switching out of CPU; will block on I/O until time " << arrival_time + t << "ms [Q"; print_queue(queue);}
        // wait_time += arrival_time;
    }
    void startIO(int t, priority_queue<Process> queue) {
        arrival_time += io_bursts[working_burst];
        if (t < 1000) {cout << "time " << t << "ms: Process "<< id << " switching out of CPU; will block on I/O until time " << arrival_time + t << "ms [Q"; print_queue(queue);}
    }
    void startIO_SRT(int t, priority_queue<Process> queue) {
        arrival_time += io_bursts[working_burst];
        cout << "time " << t << "ms: Process "<< id << " switching out of CPU; will block on I/O until time " << arrival_time + t << "ms [Q"; 
        print_queue(queue);
    }
    bool null() {return id == '.';}
    void nextCPU() {working_burst++;}
    int remainCPU() {return cpu_bursts[working_burst];}
    void wait() {wait_time++;}
    void comming() {arrival_time--;}
    void working() {cpu_bursts[working_burst]--;}
    bool finishCPU() {return remainCPU() == 0;}
    int bursts_left() {return init_number_of_bursts - working_burst - 1;}
    void set_tau(double lambda, double alpha) {
        int tau_temp = 1/lambda;
        tau.push_back(tau_temp);
        for (int i = 0; i < int(cpu_bursts.size()); i++)
        {
            tau_temp = ceil(alpha * cpu_bursts[i] + (1 - alpha) * tau_temp);
            tau.push_back(tau_temp);
        }
    }
    int current_tau() {return tau[working_burst];}
    int tau_rn() {return current_tau() - (init_cpu[working_burst] - remainCPU());}
 };

bool operator< (const Process& p1, const Process& p2)
{
    if (p1.tau[p1.working_burst] - (p1.init_cpu[p1.working_burst] - p1.cpu_bursts[p1.working_burst]) == p2.tau[p2.working_burst] - (p2.init_cpu[p2.working_burst] - p2.cpu_bursts[p2.working_burst]))
    {
        return p1.id > p2.id;
    }
    return p1.tau[p1.working_burst] - (p1.init_cpu[p1.working_burst] - p1.cpu_bursts[p1.working_burst]) > p2.tau[p2.working_burst] - (p2.init_cpu[p2.working_burst] - p2.cpu_bursts[p2.working_burst]);
}
bool operator== (const Process& p1, const Process& p2)
{
    return p1.id == p2.id  ;
}

double next_exp(double lambda, int ceiling);

void cpu(int n);

int main(int argc, char ** argv) 
// int main() 
{
    /*
    ifstream infile("commandline.txt");
    vector<string> a;
    string temp;
    while (infile >> temp)
    {
        a.push_back(temp);
    }
    int n = stoi(a[1]);
    long seed = stol(a[2]);
    double lambda = double(stol(a[3]));
    int ceiling = stoi(a[4]);
    int context_switch_time = stoi(a[5]);
    int alpha = stoi(a[6]);
    int time_slice = stoi(a[7]);
    bool rr_add = false;
    */

    if (argc < 8)
    {
        fprintf(stderr, "ERROR: Invalid arguments\n");
        return EXIT_FAILURE;
    }
    // argument 1: the number of processes to simulate. max 26
    int n = atoi(argv[1]);
    // if (n > 26 || n < 0)
    if (n < 0)
    {
        fprintf(stderr, "ERROR: Invalid arguments\n");
        return EXIT_FAILURE;
    }
    // argument 2: the seed for the random number generator
    long seed;
    sscanf(argv[2],"%lu",&seed);

    // argument 3: the lambda value to generate random value
    double lambda;
    sscanf(argv[3],"%lf",&lambda);

    // argument 4: upper bound
    int ceiling = atoi(argv[4]);

    // argument 5: context switch time
    int context_switch_time = atoi(argv[5]);
    if (context_switch_time % 2 != 0 || context_switch_time < 0)
    {
        fprintf(stderr, "ERROR: Invalid arguments\n");
        return EXIT_FAILURE;
    }

    // argument 6: constant alpha for SJF and SRT
    double alpha;
    sscanf(argv[6],"%lf",&alpha);

    // argument 7: time slice value for RR
    int time_slice = atoi(argv[7]);

    // argument 8: this is optional to see if the processes are 
    // added to the end or the beginning of the ready queue
    bool rr_add; // true means BEGINNING
    if (argc == 9)
    {
        if (string(argv[8]) == "BEGINNING") {rr_add = true;}
        else if (string(argv[8]) == "END") {rr_add = false;}
        else 
        {
            fprintf(stderr, "ERROR: Invalid arguments\n");
            return EXIT_FAILURE;
        }
    }
    else
    {
        rr_add = false;
    }

    vector<Process> p;
    srand48(seed);
    double cpu_total = 0;
    int count = 0;
    for ( int i = 0 ; i < n ; i++ )
    {
        int a_time = floor(next_exp(lambda, ceiling));
        int n_bursts = ceil(drand48() * 100);
        vector<int> temp_cpu;
        vector<int> temp_io;
        for (int j = 0; j < int(n_bursts); j++)
        {
            double cpu_burst = ceil(next_exp(lambda, ceiling));
            temp_cpu.push_back(cpu_burst);
            cpu_total += cpu_burst;
            count++;
            // printf("cpu_burst = %f\n", cpu_burst);
            if (j < n_bursts - 1)
            {
                double io_burst = ceil(next_exp(lambda, ceiling)) * 10;
                temp_io.push_back(io_burst);    
            }
            // printf("io_burst = %f\n", io_burst);
        }
        p.push_back(Process(i, a_time, n_bursts, temp_cpu, temp_io, lambda, alpha));
    }

    for ( int i = 0 ; i < int(p.size()) ; i++ )
    {
        printf("Process %c [NEW] (arrival time %d ms) %d CPU burst%s\n", p[i].id, p[i].arrival_time, p[i].init_number_of_bursts, (p[i].init_number_of_bursts == 1)?"":"s");
        // int total_c = 0;
        // int total_i = 0;
        // for (int j = 0; j < p[i].cpu_bursts.size(); j++)
        // {
        //     total_c += p[i].cpu_bursts[j];
        // }
        // for (int j = 0; j < p[i].io_bursts.size(); j++)
        // {
        //     total_i += p[i].io_bursts[j];
        // }
        // cout << total_c << ", " << total_i << endl;
    }
    
    cout << "time 0ms: Simulator started for FCFS [Q <empty>]" << endl; 
    vector<double> lst = fcfs(p, n, context_switch_time);
    printf("time %dms: Simulator ended for FCFS [Q <empty>]\n\n", int(lst[0]));
    fstream file;
    file.open ("simout.txt", fstream::out | fstream::trunc);

    file << "Algorithm FCFS\n";
    file << fixed << setprecision(3);
    file << "-- average CPU burst time: " << cpu_total / count << " ms\n";
    file << "-- average wait time: " << lst[1] << " ms\n";
    file << "-- average turnaround time: " << cpu_total / count + lst[1] + context_switch_time << " ms\n";
    file << "-- total number of context switches: " << int(lst[2]) << endl;;
    file << "-- total number of preemptions: " << int(lst[3]) << endl;;
    file << "-- CPU utilization: " << lst[4] << "%\n";

    lst = sjf(p, n, context_switch_time, lambda);
    printf("time %dms: Simulator ended for SJF [Q <empty>]\n\n", int(lst[0]));
    file << "Algorithm SJF\n";
    // file << fixed << setprecision(3);
    file << "-- average CPU burst time: " << cpu_total / count << " ms\n";
    file << "-- average wait time: " << lst[1] << " ms\n";
    file << "-- average turnaround time: " << cpu_total / count + lst[1] + context_switch_time << " ms\n";
    file << "-- total number of context switches: " << int(lst[2]) << endl;;
    file << "-- total number of preemptions: " << int(lst[3]) << endl;;
    file << "-- CPU utilization: " << lst[4] << "%\n";

    lst = srt(p, n, context_switch_time, lambda);
    printf("time %dms: Simulator ended for SRT [Q <empty>]\n\n", int(lst[0]));
    file << "Algorithm SRT\n";
    // file << fixed << setprecision(3);
    file << "-- average CPU burst time: " << cpu_total / count << " ms\n";
    file << "-- average wait time: " << lst[1] << " ms\n";
    file << "-- average turnaround time: " << lst[2] << " ms\n";
    file << "-- total number of context switches: " << int(lst[3]) << endl;;
    file << "-- total number of preemptions: " << int(lst[4]) << endl;;
    file << "-- CPU utilization: " << lst[5] << "%\n";

    for(int i = 0 ; i < int(p.size()) ; i++ )
    {
        printf("Process %c [NEW] (arrival time %d ms) %d CPU burst%s\n", p[i].id, p[i].arrival_time, p[i].init_number_of_bursts, (p[i].init_number_of_bursts == 1)?"":"s");
    }
    cout << "time 0ms: Simulator started for RR with time slice " << time_slice << "ms and rr_add to " << ((rr_add) ? "BEGINNING" : "END") << " [Q <empty>]" << endl;
    lst = rr(p, n, time_slice, rr_add, context_switch_time);
    printf("time %dms: Simulator ended for RR [Q <empty>]\n", int(lst[0]));
    file << "Algorithm RR\n";
    // file << fixed << setprecision(3);
    file << "-- average CPU burst time: " << cpu_total / count << " ms\n";
    file << "-- average wait time: " << lst[1] << " ms\n";
    file << "-- average turnaround time: " << lst[2] << " ms\n";
    file << "-- total number of context switches: " << int(lst[3]) << endl;;
    file << "-- total number of preemptions: " << int(lst[4]) << endl;;
    file << "-- CPU utilization: " << lst[5] << "%\n";

    file.close();

    return EXIT_SUCCESS;
}

double next_exp(double lambda, int ceiling)
{
    while(1)
    {
        double r = drand48();
        double x = -log( r ) / lambda;
        if ( x <= ceiling ) { return x; }
    }
}

void print_queue(vector<Process> queue)
{
    if (queue.size() == 0) {cout << " <empty>";}
    for (int i = int(queue.size()) - 1; i >= 0; i--) {
        cout << " " << queue[i].id;
    }
    cout << "]" << endl;
}

void print_queue(priority_queue<Process> queue)
{
    if (queue.size() == 0) {cout << " <empty>";}
    while (!queue.empty()) {
        cout << " " << queue.top().id;
        queue.pop();
    }
    cout << "]" << endl;
}

bool compare(const Process & l, const Process & r) //(2)
{
    return l.id < r.id;
}

int sum_up(vector<int> lst)
{
    int total = 0;
    for (int i = 0; i < int(lst.size()); i++)
    {
        total += lst[i];
    }
    return total;
}

// First-come-first-served (FCFS)
// The FCFS algorithm is a non-preemptive algorithm in which processes simply line up in the ready
// queue, waiting to use the CPU. This is your baseline algorithm (and could be implemented as RR
// with an “infinite” time slice).
vector<double> fcfs(vector<Process> processes, int n, int context_time)
{
    vector<Process> waitQueue = processes;
    vector<Process> queue;
    vector<Process> done;
    Process working_process;
    int num_terminated = 0;
    int t = 0;
    int n_preemp = 0;
    double total_working = 0;
    int total_context = 0;
    int context_wait = context_time/2;
    Process null_object;
    bool isWaiting = true;
    while(1)
    {
        sort(waitQueue.begin(), waitQueue.end(), compare);

        if (context_wait == 0 && isWaiting)
        {
            if (t < 1000) {cout << "time " << t << "ms: Process "<< working_process.id << " started using the CPU " << ((working_process.preempt) ? "with " : "for ") << working_process.remainCPU() << "ms burst" << ((working_process.preempt) ? " remaining" : "") << " [Q"; print_queue(queue);}
            isWaiting = false;
        }

        if (!working_process.null())
        {
            // if the working cpu burst finishes before timeslice
            if (working_process.finishCPU())
            {
                if (working_process.bursts_left() == 0) //all cpu bursts are over
                {
                    cout << "time " << t << "ms: Process "<< working_process.id << " terminated [Q";
                    print_queue(queue);
                    done.push_back(working_process);
                    num_terminated++;
                    if (num_terminated == n) {break;}
                }
                else // still more bursts to finish
                {
                    if (t < 1000) {cout << "time " << t << "ms: Process "<< working_process.id << " completed a CPU burst; " << working_process.bursts_left() << " burst" << ((working_process.bursts_left() == 1) ? "" : "s") << " to go [Q"; print_queue(queue); }
                    working_process.preempt = false;
                    working_process.arrival_time = context_time/2;
                    // starting I/O
                    working_process.startIO(t, queue);
                    //increment working_burst. switch to next cpu burst after I/O
                    working_process.nextCPU();
                    //push the process back to waitQueue
                    waitQueue.push_back(working_process);
                    //currently no working process here
                }
                working_process = null_object;
                context_wait = context_time;
                total_context++;
            }
        }

        // I/O burst completion
        for (int i = 0; i < int(waitQueue.size()); i++)
        {
            // cout << "checking waiting process done" << endl;
            if (waitQueue[i].arrival_time == 0 && !waitQueue[i].firstTime) // loop through waitQueue
            {
                queue.insert(queue.begin(), waitQueue[i]);
                if (t < 1000) {cout << "time " << t << "ms: Process "<< waitQueue[i].id << " completed I/O; placed on ready queue [Q"; print_queue(queue);}
                waitQueue.erase(waitQueue.begin()+i);
                i--;
            }
        }
        // New Process Arrival
        for (int i = 0; i < int(waitQueue.size()); i++)
        {
            // cout << "checking waiting process done" << endl;
            if (waitQueue[i].arrival_time == 0 && waitQueue[i].firstTime) // loop through waitQueue
            {
                waitQueue[i].firstTime = false;
                queue.insert(queue.begin(), waitQueue[i]);
                if (t < 1000) {cout << "time " << t << "ms: Process "<< waitQueue[i].id << " arrived; placed on ready queue [Q"; print_queue(queue);}
                waitQueue.erase(waitQueue.begin()+i);
                i--;
            }
        }


        if (working_process.null() && queue.size() > 0 && context_wait <= 2)
        {
            working_process = queue.back();
            queue.pop_back();
            isWaiting = true;
        }

        if (working_process.null() && !isWaiting)
        {
            if (queue.size() > 0) { context_wait--; }
            else { if (context_wait > 2) {context_wait--;} }
        }
        if (!working_process.null() && context_wait != 0) {context_wait--;}

        //increment 1ms for all processes
        if (!working_process.null() && !isWaiting)
        {
            working_process.working();
            total_working++;
            // printf("time %dms working process is %c\n", t, working_process->id);
        }

        for (int i=0; i < int(queue.size()); i++)
        {
            queue[i].wait();
        }
        for (int i=0; i < int(waitQueue.size()); i++)
        {
            waitQueue[i].comming();
        }

        t++;
    }
    t += context_time / 2;
    total_context++;

    double total_wait = 0;
    int count = 0;
    for (int i = 0; i < int(done.size()); i++)
    {
        total_wait += done[i].wait_time;
        count += done[i].init_number_of_bursts;
        // printf("%f\n", total_wait);
    }

    vector<double> items;
    items.push_back(t);
    // printf("%f\n", total_wait);
    items.push_back((total_wait - n_preemp * context_time / 2)/count);
    items.push_back(total_context);
    items.push_back(n_preemp);
    items.push_back(total_working/t * 100);
    return items;
}

// In SJF, processes are stored in the ready queue in order of priority based on their anticipated CPU
// burst times. More specifically, the process with the shortest CPU burst time will be selected as the
// next process executed by the CPU.
vector<double> sjf(vector<Process> processes, int n, int context_time, int lambda)
{
    for ( int i = 0 ; i < int(processes.size()) ; i++ )
    {
        printf("Process %c [NEW] (arrival time %d ms) %d CPU burst%s (tau %dms)\n", processes[i].id, processes[i].arrival_time, processes[i].init_number_of_bursts, (processes[i].init_number_of_bursts == 1)?"":"s", processes[i].current_tau());
    }
    cout << "time 0ms: Simulator started for SJF [Q <empty>]" << endl; 
    vector<Process> done;
    vector<Process> waitQueue = processes;
    priority_queue<Process> queue;
    Process working_process, null_object;
    int num_terminated = 0;
    int t = 0;
    int context_wait= context_time/2;
    double total_working = 0;
    int total_context = 0;
    bool isWaiting = true;
    while(1)
    {
        sort(waitQueue.begin(), waitQueue.end(), compare);

        if (context_wait == 0 && isWaiting)
        {
            if (t < 1000) {cout << "time " << t << "ms: Process "<< working_process.id << " (tau " << working_process.current_tau() << "ms) started using the CPU for " << working_process.remainCPU() << "ms burst [Q"; print_queue(queue);}
            isWaiting = false;
        }

        //CPU Burst Completion
        if (!working_process.null())
        {
            // if the working cpu burst finishes
            if (working_process.finishCPU())
            {
                if (working_process.bursts_left() == 0) //all cpu bursts are over
                {
                    cout << "time " << t << "ms: Process "<< working_process.id << " terminated [Q";
                    print_queue(queue);
                    done.push_back(working_process);
                    num_terminated++;
                    if (num_terminated == n) {break;}
                }
                else // still more bursts to finish
                {
                    if (t < 1000) {cout << "time " << t << "ms: Process "<< working_process.id << " (tau " << working_process.current_tau() << "ms) completed a CPU burst; " << working_process.bursts_left() << " burst" << ((working_process.bursts_left() == 1) ? "" : "s") << " to go [Q"; print_queue(queue); }
                    working_process.working_burst++;
                    if (t < 1000) {cout << "time " << t << "ms: Recalculated tau (" << working_process.current_tau() << "ms) for process " << working_process.id << " [Q"; print_queue(queue); }
                    working_process.working_burst--;
                    working_process.arrival_time = context_time/2;
                    // starting I/O
                    working_process.startIO(t, queue);
                    //increment working_burst. switch to next cpu burst after I/O
                    working_process.nextCPU();
                    //push the process back to waitQueue
                    waitQueue.push_back(working_process);
                    //currently no working process here
                }
                working_process = null_object;
                context_wait = context_time;
                total_context++;
            }
        }

        // I/O burst completion
        for (int i=0; i < int(waitQueue.size()); i++)
        {
            // cout << "checking waiting process done" << endl;
            if (waitQueue[i].arrival_time == 0 && !waitQueue[i].firstTime) // loop through waitQueue
            {
                queue.push(waitQueue[i]);
                if (t < 1000) {cout << "time " << t << "ms: Process "<< waitQueue[i].id << " (tau " << waitQueue[i].current_tau() << "ms) completed I/O; placed on ready queue [Q"; print_queue(queue);}
                waitQueue.erase(waitQueue.begin()+i);
                i--;
            }
        }

        // new process arrivals
        for (int i=0; i < int(waitQueue.size()); i++)
        {
            // cout << "checking waiting process done" << endl;
            if (waitQueue[i].arrival_time == 0 && waitQueue[i].firstTime) // initial arrival
            {
                waitQueue[i].firstTime = false;
                queue.push(waitQueue[i]);
                if (t < 1000) {cout << "time " << t << "ms: Process "<< waitQueue[i].id << " (tau " << waitQueue[i].current_tau() << "ms) arrived; placed on ready queue [Q"; print_queue(queue);}
                waitQueue.erase(waitQueue.begin()+i);
                i--;
            }
        }
        
        // new process starts
        if (working_process.null() && queue.size() > 0 && context_wait <= 2)
        {
            working_process = queue.top();
            queue.pop();
            isWaiting = true;
        }

        if (working_process.null() && !isWaiting)
        {
            if (queue.size() > 0) { context_wait--; }
            else { if (context_wait > 2) {context_wait--;} }
        }
        if (!working_process.null() && context_wait != 0) {context_wait--;}

        //increment 1ms for all processes
        if (!working_process.null() && !isWaiting)
        {
            working_process.working();
            total_working++;
            // printf("time %dms working process is %c\n", t, working_process->id);
        }
        priority_queue<Process> temp;
        while (!queue.empty())
        {
            Process p = queue.top();
            queue.pop();
            p.wait();
            temp.push(p);
        }
        queue = temp;
        for (int i=0; i < int(waitQueue.size()); i++)
        {
            waitQueue[i].comming();
        }
        t++;
    }
    t += context_time / 2;
    total_context++;

    double total_wait = 0;
    int count = 0;
    for (int i = 0; i < int(done.size()); i++)
    {
        total_wait += done[i].wait_time;
        count += done[i].init_number_of_bursts;
    }
    vector<double> items; 
    items.push_back(t);
    items.push_back(total_wait/count);
    items.push_back(total_context);
    items.push_back(0);
    items.push_back(total_working/t * 100);
    return items;
}

// The SRT algorithm is a preemptive version of the SJF algorithm. In SRT, when a process arrives,
// before it enters the ready queue, if it has a CPU burst time that is less than the remaining time of
// the currently running process, a preemption occurs. When such a preemption occurs, the currently
// running process is added back to the ready queue.
vector<double> srt(vector<Process> processes, int n, int context_time, int lambda)
{
    double total_burst = 0;
    for ( int i = 0 ; i < int(processes.size()) ; i++ )
    {
        total_burst += sum_up(processes[i].cpu_bursts);
        printf("Process %c [NEW] (arrival time %d ms) %d CPU burst%s (tau %dms)\n", processes[i].id, processes[i].arrival_time, processes[i].init_number_of_bursts, (processes[i].init_number_of_bursts == 1)?"":"s", processes[i].current_tau());
    }
    
    cout << "time 0ms: Simulator started for SRT [Q <empty>]" << endl;
    vector<Process> waitQueue = processes, done;
    priority_queue<Process> queue;
    Process working_process, null_object;
    int num_terminated = 0;
    int t = 0;
    int n_preemp = 0;
    double total_working = 0;
    int total_context = 0;
    int context_wait= context_time/2;
    int extra_wait = context_wait;
    bool isWaiting = true;
    bool waiting_for_statement = false;
    bool isPreempt = false;
    int a = 0;    
    while(1)
    {
        sort(waitQueue.begin(), waitQueue.end(), compare);

        if (context_wait == 0 && isWaiting)
        {
            if (t < 1000) {
                cout << "time " << t << "ms: Process "<< working_process.id << " (tau " << working_process.current_tau() << "ms) started using the CPU with " << working_process.remainCPU() << "ms burst remaining [Q"; 
                print_queue(queue);
            }
            isWaiting = false;
            if (waiting_for_statement)
            {
                working_process.preempt = true;
                Process n_working = queue.top();
                if (t < 1000) {
                    cout << "time " << t << "ms: Process "<< n_working.id << " (tau " << n_working.current_tau() << "ms) will preempt " << working_process.id << " [Q";
                    print_queue(queue);
                }
                queue.pop();
                queue.push(working_process);
                working_process = n_working;
                context_wait = context_time;
                extra_wait += context_time;
                total_context++;
                n_preemp++;
                isWaiting = true;
                waiting_for_statement = false;
            }
        }

        //CPU Burst Completion
        if (!working_process.null())
        {
            // if the working cpu burst finishes
            if (working_process.finishCPU())
            {
                if (working_process.bursts_left() == 0) //all cpu bursts are over
                {
                    cout << "time " << t << "ms: Process "<< working_process.id << " terminated [Q";
                    print_queue(queue);
                    num_terminated++;
                    done.push_back(working_process);
                    if (num_terminated == n) {break;}
                }
                else // still more bursts to finish
                {
                    if (t < 1000) {
                        cout << "time " << t << "ms: Process "<< working_process.id << " (tau " << working_process.current_tau() << "ms) completed a CPU burst; " << working_process.bursts_left() << " burst" << ((working_process.bursts_left() == 1) ? "" : "s") << " to go [Q";
                        print_queue(queue);
                    }
                    working_process.working_burst++;
                    if (t < 1000) {
                        cout << "time " << t << "ms: Recalculated tau (" << working_process.current_tau() << "ms) for process " << working_process.id << " [Q";
                        print_queue(queue);
                    }
                    working_process.working_burst--;
                    working_process.arrival_time = context_time/2;
                    // starting I/O
                    working_process.startIO(t, queue);
                    //increment working_burst. switch to next cpu burst after I/O
                    working_process.nextCPU();
                    //push the process back to waitQueue
                    waitQueue.push_back(working_process);
                    //currently no working process here
                }
                working_process = null_object;
                context_wait = context_time;
                total_context++;
                extra_wait += context_time;
            }
        }

        // I/O burst completion
        for (int i=0; i < int(waitQueue.size()); i++)
        {
            // cout << "checking waiting process done" << endl;
            if (waitQueue[i].arrival_time == 0 && !waitQueue[i].firstTime) // loop through waitQueue
            {
                if (t < 1000) {
                    cout << "time " << t << "ms: Process "<< waitQueue[i].id << " (tau " << waitQueue[i].current_tau() << "ms) completed I/O; ";
                }
                if (!working_process.null() && waitQueue[i].current_tau() < working_process.tau_rn() && !isWaiting)
                {
                    queue.push(waitQueue[i]);
                    working_process.preempt = true;
                    if (t < 1000) {
                        cout << "preempting " << working_process.id << " [Q";
                        print_queue(queue);
                    }
                    isPreempt = true;
                    queue.push(working_process);
                    working_process = waitQueue[i];
                    context_wait = context_time;
                    extra_wait += context_time;
                    total_context++;
                    n_preemp++;
                    isWaiting = true;
                }
                else
                {
                    if (!working_process.null() && waitQueue[i].current_tau() < working_process.tau_rn()) {waiting_for_statement = true;}
                    queue.push(waitQueue[i]);
                    if (t < 1000) {
                        cout << "placed on ready queue [Q";
                        print_queue(queue);
                    }   
                }
                waitQueue.erase(waitQueue.begin()+i);
                i--;
            }
        }

        // new process arrivals
        for (int i=0; i < int(waitQueue.size()); i++)
        {
            if (waitQueue[i].arrival_time == 0 && waitQueue[i].firstTime) // initial arrival
            {
                if (t < 1000) {
                    cout << "time " << t << "ms: Process "<< waitQueue[i].id << " (tau " << waitQueue[i].current_tau() << "ms) arrived; ";
                }
                waitQueue[i].firstTime = false;
                if (!working_process.null())
                {
                    if (waitQueue[i].current_tau() < working_process.tau_rn() && !isWaiting)
                    {
                        working_process.preempt = true;
                        queue.push(waitQueue[i]);
                        if (t < 1000) {
                            cout << "preempting " << working_process.id << " [Q";
                            print_queue(queue);
                        }
                        isPreempt = true;
                        queue.push(working_process);
                        working_process = waitQueue[i];
                        context_wait = context_time;
                        extra_wait += context_time;
                        total_context++;
                        n_preemp++;
                        isWaiting = true;
                    }
                    else
                    {
                        if (waitQueue[i].current_tau() < working_process.tau_rn()) {waiting_for_statement = true;}
                        queue.push(waitQueue[i]);
                        if (t < 1000) {
                            cout << "placed on ready queue [Q";
                            print_queue(queue);
                        }
                    }
                }
                else
                {
                    queue.push(waitQueue[i]);
                    if (t < 1000) {
                        cout << "placed on ready queue [Q";
                        print_queue(queue);
                    }
                }
                waitQueue.erase(waitQueue.begin()+i);
                i--;
            }
        }
        
        // new process starts
        if (working_process.null() && queue.size() > 0 && context_wait <= 2)
        {
            working_process = queue.top();
            queue.pop();
            isWaiting = true;
        }

        if (working_process.null() && !isWaiting)
        {
            if (queue.size() > 0) { context_wait--; }
            else { if (context_wait > 2) {context_wait--;} }
        }
        if (!working_process.null() && context_wait != 0) {context_wait--;}
        if (isPreempt && context_wait < 2) 
        {
            isPreempt = false; 
            if (!(working_process == queue.top())) 
            {working_process = queue.top(); waiting_for_statement = false;}
            queue.pop();
            a++;
        }

        //increment 1ms for all processes
        if (!working_process.null() && !isWaiting)
        {
            working_process.working();
            total_working++;
            // working_process.tau_working();
            // printf("time %dms working process is %c\n", t, working_process->id);
        }
        priority_queue<Process> temp;
        while (!queue.empty())
        {
            Process p = queue.top();
            queue.pop();
            p.wait();
            temp.push(p);
        }
        queue = temp;
        for (int i=0; i < int(waitQueue.size()); i++)
        {
            waitQueue[i].comming();
        }

        t++;
    }
    t += context_time / 2;
    extra_wait += context_time / 2;
    total_context++;

    double total_wait = 0;
    int count = 0;
    for (int i = 0; i < int(done.size()); i++)
    {
        total_wait += done[i].wait_time;
        count += done[i].init_number_of_bursts;
        // printf("%f\n", total_wait);
    }

    vector<double> items;
    items.push_back(t);
    items.push_back((total_wait - (n_preemp - a) * context_time / 2)/count); //wait time
    items.push_back((total_wait - (n_preemp - a) * context_time/2 + total_burst + extra_wait)/count); //turnaround
    items.push_back(total_context); //number of context switch
    items.push_back(n_preemp); // number of preemptions
    items.push_back(total_working/t * 100); // utilization
    return items;
}

// The RR algorithm is essentially the FCFS algorithm with predefined time slice tslice. Each process
// is given tslice amount of time to complete its CPU burst. If this time slice expires, the process is
// preempted and added to the end of the ready queue (though see the rradd parameter described below).
// If a process completes its CPU burst before a time slice expiration, the next process on the ready
// queue is immediately context-switched in to use the CPU.
vector<double> rr(vector<Process> processes, int n, int timeslice, bool rr_add, int context_time)
{
    double total_burst = 0;
    for (int i = 0; i < int(processes.size()); i++)
    {
        total_burst += sum_up(processes[i].cpu_bursts);
    }
    vector<Process> waitQueue = processes;
    vector<Process> queue;
    vector<Process> done;
    Process working_process;
    int num_terminated = 0;
    int time_spent = 0;
    int t = 0;
    int n_preemp = 0;
    double total_working = 0;
    int total_context = 0;
    int context_wait = context_time/2;
    int extra_wait = context_wait;
    Process null_object;
    bool isWaiting = true;
    while(1)
    {
        sort(waitQueue.begin(), waitQueue.end(), compare);

        if (context_wait == 0 && isWaiting)
        {
            if (t < 1000) {cout << "time " << t << "ms: Process "<< working_process.id << " started using the CPU " << ((working_process.preempt) ? "with " : "for ") << working_process.remainCPU() << "ms burst" << ((working_process.preempt) ? " remaining" : "") << " [Q";print_queue(queue);}
            isWaiting = false;
            time_spent = 0;
        }

        if (!working_process.null())
        {
            // if the working cpu burst finishes before timeslice
            if (working_process.finishCPU())
            {
                if (working_process.bursts_left() == 0) //all cpu bursts are over
                {
                    cout << "time " << t << "ms: Process "<< working_process.id << " terminated [Q";
                    print_queue(queue);
                    done.push_back(working_process);
                    num_terminated++;
                    if (num_terminated == n) {break;}
                }
                else // still more bursts to finish
                {
                    if (t < 1000) {cout << "time " << t << "ms: Process "<< working_process.id << " completed a CPU burst; " << working_process.bursts_left() << " burst" << ((working_process.bursts_left() == 1) ? "" : "s") << " to go [Q";print_queue(queue);}
                    working_process.preempt = false;
                    working_process.arrival_time = context_time/2;
                    // starting I/O
                    working_process.startIO(t, queue);
                    //increment working_burst. switch to next cpu burst after I/O
                    working_process.nextCPU();
                    //push the process back to waitQueue
                    waitQueue.push_back(working_process);
                    //currently no working process here
                }
                working_process = null_object;
                context_wait = context_time;
                total_context++;
                extra_wait += context_time;
            }

            // check if the process is expired
            else if (time_spent == timeslice && context_wait == 0)
            {
                if (queue.size() == 0)
                {
                    if (t < 1000) {cout << "time " << t << "ms: Time slice expired; no preemption because ready queue is empty [Q <empty>]" << endl;}
                }
                else
                {
                    if (t < 1000) {cout << "time " << t << "ms: Time slice expired; process "<< working_process.id << " preempted with " << working_process.remainCPU() << "ms to go [Q"; print_queue(queue);}
                    working_process.preempt = true;
                    queue.insert(queue.begin(), working_process);
                    //currently no working process here
                    working_process = null_object;
                    context_wait = context_time;
                    extra_wait += context_time;
                    // total_context += context_time / 2;
                    total_context++;
                    n_preemp++;
                }
                time_spent = 0;
            }
        }

        // I/O burst completion
        for (int i=0; i < int(waitQueue.size()); i++)
        {
            // cout << "checking waiting process done" << endl;
            if (waitQueue[i].arrival_time == 0 && !waitQueue[i].firstTime) // loop through waitQueue
            {
                if (rr_add) { queue.push_back(waitQueue[i]); }
                else { queue.insert(queue.begin(), waitQueue[i]); }
                if (t < 1000) {cout << "time " << t << "ms: Process "<< waitQueue[i].id << " completed I/O; placed on ready queue [Q"; print_queue(queue);}
                waitQueue.erase(waitQueue.begin()+i);
                i--;
            }
        }
        // New Process Arrival
        for (int i=0; i < int(waitQueue.size()); i++)
        {
            // cout << "checking waiting process done" << endl;
            if (waitQueue[i].arrival_time == 0 && waitQueue[i].firstTime) // loop through waitQueue
            {
                waitQueue[i].firstTime = false;
                if (rr_add) { queue.push_back(waitQueue[i]); }
                else { queue.insert(queue.begin(), waitQueue[i]); }
                if (t < 1000) {cout << "time " << t << "ms: Process "<< waitQueue[i].id << " arrived; placed on ready queue [Q"; print_queue(queue);}
                waitQueue.erase(waitQueue.begin()+i);
                i--;
            }
        }


        if (working_process.null() && queue.size() > 0 && context_wait <= 2)
        {
            working_process = queue.back();
            queue.pop_back();
            isWaiting = true;
        }

        if (working_process.null() && !isWaiting)
        {
            if (queue.size() > 0) { context_wait--; }
            else { if (context_wait > 2) {context_wait--;} }
        }
        if (!working_process.null() && context_wait != 0) {context_wait--;}

        //increment 1ms for all processes
        if (!working_process.null() && !isWaiting)
        {
            working_process.working();
            total_working++;
            // printf("time %dms working process is %c\n", t, working_process->id);
        }

        for (int i=0; i < int(queue.size()); i++)
        {
            queue[i].wait();
        }
        for (int i=0; i < int(waitQueue.size()); i++)
        {
            waitQueue[i].comming();
        }

        t++;
        time_spent++;
    }
    t += context_time / 2;
    extra_wait += context_time / 2;
    total_context++;

    double total_wait = 0;
    int count = 0;
    for (int i = 0; i < int(done.size()); i++)
    {
        total_wait += done[i].wait_time;
        count += done[i].init_number_of_bursts;
        // printf("%f\n", total_wait);
    }

    vector<double> items;
    items.push_back(t);
    // printf("%f\n", total_wait);
    items.push_back((total_wait - n_preemp * context_time / 2)/count); //wait time
    items.push_back((total_wait - n_preemp * context_time/2 + total_burst + extra_wait)/count); //turnaround
    items.push_back(total_context); //number of context switch
    items.push_back(n_preemp); // number of preemptions
    items.push_back(total_working/t * 100); // utilization
    return items;
}