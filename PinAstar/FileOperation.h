#pragma once

#ifndef _FILEOPERATION_
#define _FILEOPERATION_

#include"DefineClass.h"
#include <string>
#include <windows.h>
#include <fstream>
void ClearFile(string name);
void WriteFile(string name, double input);
void WriteFile(string name,string input);
void WriteFile_Parameters(string Name,CODE LinearBlockCode, OPER_PARA key_in_data);
void WriteFile(string name, vector<size_t> input);
//void WriteFile2(string name, vector<size_t> input);
void Show_Current_Time();

#endif
//File
