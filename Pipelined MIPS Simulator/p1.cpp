/* p1.cpp */
/* NAME: Yudai Teruyama */
#include <cassert>
#include <fstream> //For file I/O
#include <iostream> //For std::cerr
#include <string> //For parsing
#include <ctime> //For std::time_t
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include "Instruction.h"
using namespace std;

#define EXIT_SUCCESS 0 
#define EXIT_FAILURE 1 

void NonForward(vector <string> registers, vector <int> dataHazard, unsigned int size, vector <Instruction> file, vector<int> loop_pos, vector<string> loop_name);
void Forward(vector <string> registers, vector <int> dataHazard, unsigned int size, vector <Instruction> file, vector<int> loop_pos, vector<string> loop_name);
void move_on(vector <Instruction> &file, vector <string> &registers);
void Branch(unsigned int cycle, unsigned int &i, unsigned int &instNum, unsigned int &size, vector <Instruction> &file, vector<int> loop_pos, vector<string> loop_name);

void startCycle() {
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "CPU Cycles ===>     1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16" << endl;
}

void endCycle() {
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "END OF SIMULATION" << endl;
}

void print_vals(vector <string> registers) { //function to print the info of all variables
	int col = 0;
	for (unsigned int i=0; i<18; i++) {
		if (i < 8) {
			cout << "$s" << i << " = ";
		}
		else {
			cout << "$t" << i - 8 << " = ";
		}
		if (col == 3 || i==17) {
			cout << registers[i] << endl;
			col = 0;
		}
		else {
			cout.width(14); 
			cout << registers[i];
			col++;
		}
	}
}

void hazardsCheck(vector <Instruction> file, vector <int> &dataHazard) {
	//check if there will be any possible datahazard 
	//make sure to put how many stages the program needs to wait
	vector<string> prev{" ", " ", " ", " ", " ", " ", " ", " ", " ", " "};
	int num = 0;
	unsigned int i, j, k;
	for (i=0; i<file.size(); i++) {
		prev[num] = file[i].getr1();
		num++;
		if (i>=3) {j=i-2;}
		else {j=0;}
		for (k=j; k<i; k++) {
			if ((prev[k] == file[i].getr2()) || (prev[k] == file[i].getr3())) {
				dataHazard[i] = 4 + k - (i + 2);
			}
		}
	}
}

int main(int argc, char** argv){
	if(argc != 3){
		cerr << "Correct usage is " << argv[0] << " [input file] [output file]" << endl;
		return -1;
	}

	string mode = argv[1];
	//open a file
	ifstream infile(argv[2]);
	if(!infile){
		cerr << "Could not open " << argv[1] << " for reading!";
		return -1;
	}
	string inst, RsRtRd;
	vector <Instruction> file;
	vector<int> loop_pos;
	vector<string> loop_name;
	vector <string> registers{"0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0"}; 
	vector <int> dataHazard = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	unsigned int size = 0, loop_start = 0;
	Instruction LINE;
	//read file and create class for each instruction
	while(infile >> inst){
		if (inst[inst.length()-1] == ':') {
			inst.erase(inst.end()-1);
			loop_name.push_back(inst);
			if (loop_start != 0) {
				loop_pos.push_back(loop_start);
			}
			loop_start = size;
		}
		else {
			infile >> RsRtRd;
			LINE = Instruction(inst, RsRtRd);
			file.push_back(LINE);
			size++;
		}
	}
	if (loop_start != 0) {
		loop_pos.push_back(loop_start);
	}

	hazardsCheck(file, dataHazard);
	//the simulation starts
	cout << "START OF SIMULATION ";
	if (mode == "N") {NonForward(registers, dataHazard, size, file, loop_pos, loop_name);}
	if (mode == "F") {Forward(registers, dataHazard, size, file, loop_pos, loop_name);}
	endCycle();
	return EXIT_SUCCESS;
}

int move_on(Instruction &inst, vector <string> &registers, unsigned int cycle, unsigned int instNum, int nopnop) {
	//move each instruction by one stage
	vector<string> temp = inst.getCycle();
	unsigned int position = cycle-1;
	if (temp[position] != ".") {;}
	else if (inst.getProcess() == 0) {temp[position] = "IF";}
	else if (inst.getNop() > 0) {temp[position] = temp[position-1]; inst.minusNop();}
	else if (temp[position-1] == "IF") {temp[position] = "ID";}
	else if (temp[position-1] == "ID") {temp[position] = "EX";}
	else if (temp[position-1] == "EX") {temp[position] = "MEM";}
	else if (temp[position-1] == "MEM") {
		temp[position] = "WB";
		if (inst.getOP() == "bne" || inst.getOP() == "beq") {
			inst.copy_cycle(temp);
			inst.moveOn();
			return 1;
		}
		else {
			executeInstruction(inst, registers);
		}
	}
	else if (temp[position-1] == "*" && inst.getStar() < 2) {temp[position] = "*"; inst.starPlus();}
	inst.copy_cycle(temp);
	inst.moveOn();
	return 0;
}

void NonForward(vector <string> registers, vector <int> dataHazard, unsigned int size, vector <Instruction> file, vector<int> loop_pos, vector<string> loop_name) {
	//function to execute instructions without forwarding
	cout << "(no forwarding)" << endl;
	unsigned int cycle = 1, instNum = 1, i;
	int nop_count, nopnop = 0, branch_check;
	while (cycle < 17 && cycle < size + 5) {
		nop_count = 0;
		startCycle();
		for (i=0; i<instNum; i++) {
			if (file[i].getProcess() >= 2 && dataHazard[i] != -1) {
				if (dataHazard[i] != -1 && nop_count == 0) {// insert nop class if needed
					Instruction nop;
					vector<Instruction>::iterator newIt;
					vector<int>::iterator newIt1;
					for (unsigned int j=i; j<instNum; j++) {
						file[j].plusNop();
					}
					nop = Instruction("nop", " ");
					nop.copy_cycle(file[i].getCycle());
					nop.stage_change("*", cycle-1);
					nop.moveOn();
					newIt = file.insert(file.begin() + i, nop);
					instNum++;
					nop_count++;
					nopnop++;
					
					dataHazard[i] -= 1;
					newIt1 = dataHazard.insert(dataHazard.begin() + i, -1);

					cout << file[i];
					i++;
					branch_check = move_on(file[i], registers, cycle, i, nopnop);
				}
				cout << file[i];
			}
			else {
				//check if there is any control hazard
				branch_check = move_on(file[i], registers, cycle, i, nopnop);
				cout << file[i];
				if (branch_check == 1) {
				if(branch(file[i], registers)) {
					Branch(cycle, i, instNum, size, file, loop_pos, loop_name);
				}
			}
			}
		}
		cout << endl;
		print_vals(registers);
		cycle++;
		if (nop_count) {size+=nop_count;}
		if (instNum < size) {instNum++;}
	}
}
void Forward(vector <string> registers, vector <int> dataHazard, unsigned int size, vector <Instruction> file, vector<int> loop_pos, vector<string> loop_name) {
	//function to execute instructions with forwarding
	cout << "(forwarding)" << endl;
	unsigned int cycle = 1, instNum = 1, i;
	int branch_check;
	vector<string> temp;
	string op, loop;
	vector<string> stages;
	while (cycle < 17 && cycle < size + 5) {
		startCycle();
		for (i=0; i<instNum; i++) {
			branch_check = move_on(file[i], registers, cycle, i, 0);
			cout << file[i];
			if (branch_check == 1) {
				//no need to worry abt data hazard just control hazard
				if(branch(file[i], registers)) {
					Branch(cycle, i, instNum, size, file, loop_pos, loop_name);
				}
			}
		}
		cout << endl;
		print_vals(registers);
		cycle++;
		if (instNum < size) {instNum++;}
	}
}

void Branch(unsigned int cycle, unsigned int &i, unsigned int &instNum, unsigned int &size, vector <Instruction> &file, vector<int> loop_pos, vector<string> loop_name) {
	//if it branches...
	vector<string> temp;
	string loop;
	unsigned int j, h, start, end; 
	unsigned int num = 1;
	for (j=i+1; j<instNum; j++) {
		temp = file[j].getCycle();
		temp[cycle-1] = "*";
		file[j].copy_cycle(temp);
		file[j].moveOn();
		file[j].changeStar(num);
		num--;
		cout << file[j];
	}
	loop = file[i].getr3();
	for (h=0; h<loop_name.size(); h++) {
		if (loop == loop_name[h]) {
			start = loop_pos[h];
			break;
		}
	}
	if (h != loop_pos.size()-1) {end = loop_pos[h+1];}
	else {end = size;}
	Instruction tem;
	vector<Instruction>::iterator newIt;
	for (unsigned int x=start; x<end; x++) {
		tem = file[x];
		tem.reset();
		newIt = file.insert(file.begin()+j, tem);
		size++;
		j++;
	}
	i = instNum-1;
	instNum++;
}