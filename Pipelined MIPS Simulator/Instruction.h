#ifndef __INSTRUCTION_
#define __INSTRUCTION_
#include <string>
#include <ctime>    //Needed for std::time_t
#include <iostream>    //Needed for std::ostream
#include <vector>
#include <array>
#include <iomanip>
using namespace std;

class Instruction{
public:
    Instruction(string opcode, string registers);
    Instruction();
    Instruction(const Instruction& inst);
    void copy_cycle(vector<string> temp);
    string getOP() const { return op; }
    string getr1() const { return r1; }
    string getr2() const { return r2; }
    string getr3() const { return r3; }
    string getLine() const { return line; }
    int getProcess() const { return process; }
    int getStar() const { return star; }
    string cuurentStage() const { return cycle[process]; }
    vector<string> getCycle() const { return cycle; }
    void stage_change(string str, int pos);
    void moveOn() {process++;}
    void starPlus() {star++;}
    void changeStar(int num) {star = num;}
    void reset();
    int getNop() const {return nopcount;}
    void plusNop() {nopcount++;}
    void minusNop() {nopcount--;}


private:
    string op;
    string r1;
    string r2;
    string r3;
    string line;
    vector<string> cycle;
    int process;
    int star;
    int nopcount;
};

//Promise to declare stream output operator for a instruction 
ostream& operator<<(ostream& out, const Instruction inst);
//execute each instruction
void executeInstruction(const Instruction inst, vector <string>& registers);
//check if it branches
int branch(const Instruction inst, vector <string>& registers);
//calculte the length of each register
int range(string str, int start);

#endif