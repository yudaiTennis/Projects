#include "Instruction.h"
using namespace std;
#include <string>


Instruction::Instruction() {//initialization 
    op = " ";
    r1 = r2 = r3 = line = " ";
    process = star = nopcount = 0;
    for (unsigned int i=0; i<16; i++) {
        cycle.push_back(".");
    }
}
//Constructor for the Instruction class
Instruction::Instruction(string opcode, string registers) {
    op = opcode;
    int num = 0;
    process = star = nopcount = 0;
    line = opcode + ' ' + registers;
    for (unsigned int i=0; i<registers.length(); i++) {
        if (num == 0 && registers[i] == '$') {
            int rg = range(registers, i+1);
            r1 = registers.substr(i+1, rg);
            i += 2;
            num++;
        }
        else if (registers[i] == ',' && num == 1) {
            if (registers[i+1] == '$') {i++;}
            int rg = range(registers, i+1);
            r2 = registers.substr(i+1, rg);
            num++;
            i += rg;
        }
        else if (registers[i] == ',' && num == 2) {
            if (registers[i+1] == '$') {i++;}
            int rg = range(registers, i+1);
            r3 = registers.substr(i+1, rg);
            num++;
            i += rg;
        }
    }
    for (unsigned int i=0; i<16; i++) {
        cycle.push_back(".");
    }
}


int range(string str, int start) {
    int i=0;
    while (start != int(str.length())) {
        if (str[start] == ',') {break;}
        i++;
        start++;
    }
    return i;
}
//Prints a instruction, given an output stream
ostream& operator<<(ostream& out, const Instruction inst){
    out.fill(' ');
    out.width(20);
    out << left << inst.getLine();
    vector<string> cy = inst.getCycle();
    for (unsigned int i=0; i<cy.size(); i++) {
        if (i == 15) {out << cy[i];}
        else {out << setw(4) << cy[i];}
    }
    out << endl;    
    return out;
}
//get a value for a register
int findVal(string rg, vector <string> registers) {
    int val = 0;
    if (rg[0] == 's'){
        val = std::stoi(registers[rg[1]-'0']);
    }
    else if (rg[0] == 't'){
        val = std::stoi(registers[rg[1]-'0'+8]);
    }
    else if (rg == "zero"){
        val = 0;
    }
    else {
        val = std::stoi(rg);
    }
    return val;
}

void executeInstruction(const Instruction inst, vector <string>& registers) {
    string op = inst.getOP();
    string r1 = inst.getr1();
    string r2 = inst.getr2();
    string r3 = inst.getr3();
    int val_r1 = 0, val_r2 = findVal(r2, registers), val_r3 = findVal(r3, registers);

    if(op == "add" || op == "addi") {
        val_r1 = val_r2 + val_r3;
    }
    else if(op == "and" || op == "andi") {
        val_r1 = val_r2 & val_r3;
    }
    else if(op == "or" || op == "ori") {
        val_r1 = val_r2 | val_r3;
    }
    else if(op == "slt" || op == "slti") {
        if (val_r2 < val_r3) {val_r1 = 1;}
    }

    if (r1[0] == 's'){
        registers[r1[1]-'0'] = to_string(val_r1);
    }
    else if (r1[0] == 't'){
        registers[r1[1]-'0'+8] = to_string(val_r1);
    } 
}

int branch(const Instruction inst, vector <string>& registers) {
    string op = inst.getOP();
    string r1 = inst.getr1();
    string r2 = inst.getr2();
    string r3 = inst.getr3();
    int val_r1 = findVal(r1, registers), val_r2 = findVal(r2, registers);
    if (op == "bne") {
        if (val_r1 == val_r2) {return 0;}
        return 1;
    }
    else if (op == "beq") {
        if (val_r1 != val_r2) {return 0;}
        return 1;
    }
    return 1;
}
//copy instruction
Instruction::Instruction(const Instruction& inst){
    op = inst.op; 
    r1 = inst.r1;
    r2 = inst.r2;
    r3 = inst.r3;
    line = inst.line;
    cycle = inst.cycle;
    process = inst.process;
    star = inst.star;
    nopcount = inst.nopcount;
}

void Instruction::copy_cycle(vector<string> temp){
    cycle = temp;
}

void Instruction::stage_change(string str, int pos) {
    cycle[pos] = str;
}
//reset a cycle
void Instruction::reset() {
    vector<string> temp;
    for (unsigned int i=0; i<16; i++) {
        temp.push_back(".");
    }
    cycle = temp;
    process = star = 0;
}