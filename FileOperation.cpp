#include "PinAstar/FileOperation.h"

void Show_Current_Time()
{
	SYSTEMTIME sys;
	GetLocalTime(&sys);
	cout
		<< "\n "
		<< sys.wYear << " / "
		<< sys.wMonth << " / "
		<< sys.wDay << " , "
		<< sys.wHour << " : "
		<< sys.wMinute << " : "
		<< sys.wSecond;
}

void WriteFile_Parameters(string Name, CODE LinearBlockCode, OPER_PARA key_in_data)
{
	string name = Name;
	ClearFile(name);
	using namespace std;
	SYSTEMTIME sys;	//windows API
	GetLocalTime(&sys);

	string filename(name);
	fstream fp;
	fp.open(filename, ios::out | ios::app);
	if (!fp)
		cout << "Fail to open file: " << filename << endl;

	switch (key_in_data.Decoder_Version)
	{
	case A_Star_Original:
		fp << " \n A* I (Original)";
		break;
	case A_Star_Parity:
		fp << " \n A* Parity";
		break;
	case A_Star_Parity_f:
		fp << " \n A* Parity-f, Alpha = " << key_in_data.Alpha;
		break;

		// RSS
	case A_Star_PC:
		fp
			<< " \n A* PC-"
			<< key_in_data.Constraint_i << " ";
		break;
	case A_Star_PC_out:
		fp
			<< " \n A* PC-out-"
			<< key_in_data.Constraint_i;
		break;
	case A_Star_2_Stack:
		fp << " \n A* 2-Stack";
		break;
	case A_Star_3_Stack:
		fp << " \n A* 3-Stack";
		break;
	case A_Star_4_Stack:
		fp << " \n A* 4-Stack";
		break;

		// RSS + Parity Heuristic Function
	case A_Star_Parity_PC:
		fp
			<< " \n A* Parity + PC-"
			<< key_in_data.Constraint_i << " ";
		break;
	case A_Star_Parity_PC_out:
		fp
			<< " \n A* Parity + PC-out-"
			<< key_in_data.Constraint_i << " ";
		break;
	case A_Star_2_Parity_Stack:
		fp << " \n A* Parity + 2-Stack";
		break;
	case A_Star_3_Parity_Stack:
		fp << " \n A* Parity + 3-Stack";
		break;
	case A_Star_4_Parity_Stack:
		fp << " \n A* Parity + 4-Stack";
		break;

		// DM-I 
	case A_Star_PC_out_DM_I:
		fp
			<< " \n A* PC-out-" << key_in_data.Constraint_i
			<< "-" << key_in_data.Constraint_j << " DM-I "
			<< ", Control Level = " << key_in_data.Control_Level;
		break;
	case A_Star_2_Stack_DM_I:
		fp
			<< " \n A* 2-Stack-" << key_in_data.Constraint_j << " DM-I "
			<< ", Control Level = " << key_in_data.Control_Level;
		break;
	case A_Star_3_Stack_DM_I:
		fp
			<< " \n A* 3-Stack-" << key_in_data.Constraint_j << " DM-I "
			<< ", Control Level = " << key_in_data.Control_Level;
		break;
	case A_Star_4_Stack_DM_I:
		fp
			<< " \n A* 4-Stack-" << key_in_data.Constraint_j << " DM-I "
			<< ", Control Level = " << key_in_data.Control_Level;
		break;

		// DM-I + Parity Heuristic Function 
	case A_Star_PC_out_DM_I_Parity:
		fp
			<< " \n A* PC-out-" << key_in_data.Constraint_i
			<< "-" << key_in_data.Constraint_j
			<< " DM-I + Parity"
			<< ", Control Level = " << key_in_data.Control_Level;
		break;
	case A_Star_2_Stack_DM_I_Parity:
		fp
			<< " \n A* 2-Stack-" << key_in_data.Constraint_j
			<< " DM-I + Parity"
			<< ", Control Level = " << key_in_data.Control_Level;
		break;
	case A_Star_3_Stack_DM_I_Parity:
		fp
			<< " \n A* 3-Stack-" << key_in_data.Constraint_j
			<< " DM-I + Parity"
			<< ", Control Level = " << key_in_data.Control_Level;
		break;
	case A_Star_4_Stack_DM_I_Parity:
		fp
			<< " \n A* 4-Stack-" << key_in_data.Constraint_j
			<< " DM-I + Parity"
			<< ", Control Level = " << key_in_data.Control_Level;
		break;

		// DM-II
	case A_Star_PC_DM_II:
		fp
			<< " \n  A* PC-DM-II-"
			<< key_in_data.Constraint_i << "-" << key_in_data.Constraint_j << " "
			<< "\n Check Level = " << key_in_data.Check_Level
			<< "\n Control Level = " << key_in_data.Control_Level
			<< "\n DM Stack Size = " << key_in_data.DM_StackSize;
		break;
	case A_Star_PCout_DM_II:
		fp
			<< " \n A* PC-out-DM-II-"
			<< key_in_data.Constraint_i << "-" << key_in_data.Constraint_j << " "
			<< "\n Check Level = " << key_in_data.Check_Level
			<< "\n Control Level = " << key_in_data.Control_Level
			<< "\n DM Stack Size = " << key_in_data.DM_StackSize;
		break;
	case A_Star_2_Stack_DM_II:
		fp
			<< " \n A* 2-Stack-"
			<< key_in_data.Constraint_j << "-DM-II "
			<< "\n Check Level = " << key_in_data.Check_Level
			<< "\n Control Level = " << key_in_data.Control_Level
			<< "\n DM Stack Size = " << key_in_data.DM_StackSize;
		break;
	case A_Star_3_Stack_DM_II:
		fp
			<< " \n A* 3-Stack-"
			<< key_in_data.Constraint_j << "-DM-II "
			<< "\n Check Level = " << key_in_data.Check_Level
			<< "\n Control Level = " << key_in_data.Control_Level
			<< "\n DM Stack Size = " << key_in_data.DM_StackSize;
		break;
	case A_Star_4_Stack_DM_II:
		fp
			<< " \n A* 4-Stack-"
			<< key_in_data.Constraint_j << "-DM-II "
			<< "\n Check Level = " << key_in_data.Check_Level
			<< "\n Control Level = " << key_in_data.Control_Level
			<< "\n DM Stack Size = " << key_in_data.DM_StackSize;
		break;
	case A_Star_BMA:
		fp
			<< " \n A* BMA-"
			<< key_in_data.Constraint_j << " "
			<< "\n Check Level = " << key_in_data.Check_Level
			<< "\n Control Level = " << key_in_data.Control_Level
			<< "\n DM Stack Size = " << key_in_data.DM_StackSize;
		break;

		// OSC
	case A_Star_PC_OSC:
		fp
			<< "\n A* PC-" << key_in_data.Constraint_i
			<< "\n OSC-" << key_in_data.OSC_Alpha;
		break;
	case A_Star_PCout_OSC:
		fp
			<< "\n A* PC-out-" << key_in_data.Constraint_i
			<< "\n OSC-" << key_in_data.OSC_Alpha;
		break;
	case A_Star_2_Stack_OSC:
		fp
			<< "\n A* 2-Stack"
			<< "\n OSC-" << key_in_data.OSC_Alpha;
		break;
	case A_Star_3_Stack_OSC:
		fp
			<< "\n A* 3-Stack"
			<< "\n OSC-" << key_in_data.OSC_Alpha;
		break;
	case A_Star_4_Stack_OSC:
		fp
			<< "\n A* 4-Stack"
			<< "\n OSC-" << key_in_data.OSC_Alpha;
		break;

		// CBC
	case A_Star_PCout_CBC:
		fp
			<< "\n A* PC-out-" << key_in_data.Constraint_i
			<< "\n CBC-" << key_in_data.CBC_FlippingBit;
		break;
	case A_Star_2_Stack_CBC:
		fp
			<< "\n A* 2-Stack"
			<< "\n CBC-" << key_in_data.CBC_FlippingBit;
		break;
	case A_Star_3_Stack_CBC:
		fp
			<< "\n A* 3-Stack"
			<< "\n CBC-" << key_in_data.CBC_FlippingBit;
		break;

		// CBC + OSC
	case A_Star_PCout_CBC_OSC:
		fp
			<< "\n A* PC-out-" << key_in_data.Constraint_i
			<< "\n OSC-" << key_in_data.OSC_Alpha
			<< "\n CBC-" << key_in_data.CBC_FlippingBit;
		break;
	case A_Star_2_Stack_CBC_OSC:
		fp
			<< "\n A* 2-Stack" << key_in_data.Constraint_i
			<< "\n OSC-" << key_in_data.OSC_Alpha
			<< "\n CBC-" << key_in_data.CBC_FlippingBit;
		break;
	case A_Star_3_Stack_CBC_OSC:
		fp
			<< "\n A* 3-Stack" << key_in_data.Constraint_i
			<< "\n OSC-" << key_in_data.OSC_Alpha
			<< "\n CBC-" << key_in_data.CBC_FlippingBit;
		break;
	}


	fp
		<< " Decoding Algorithm Simulation for " << LinearBlockCode.Title << "\n"
		<< "\n -----------------------------------------------------\n\n"
		<< " Start Time : " << sys.wYear << " / " << sys.wMonth << " / " << sys.wDay << " , " << sys.wHour << " : " << sys.wMinute << "\n\n"
		<< " The Number of Blocks =  " << key_in_data.block_number << "\n"
		<< " The Number of Symbols =  " << LinearBlockCode.Row_Number *  key_in_data.block_number << "\n"
		<< " Start SNR [dB] =  " << key_in_data.SNR_dB_start << " [dB] ¡÷ Es/N0 [dB] = " << key_in_data.SNR_dB_start + 10 * log10(LinearBlockCode.Code_Rate) << " [dB]\n"
		<< " End SNR [dB] =  " << key_in_data.SNR_dB_end << " [dB] ¡÷ Es/N0 [dB] = " << key_in_data.SNR_dB_end + 10 * log10(LinearBlockCode.Code_Rate) << " [dB]\n"
		<< " Step SNR [dB] =  " << key_in_data.SNR_dB_step << " [dB]" << "\n"
		<< " Block Step = " << key_in_data.BlockStep << "\n"
		<< " Stack Size = " << key_in_data.StackSize << "\n"
		<< " Error Blocks Threshold = " << key_in_data.ErrorBlock_Thr << "\n"
		<< "\n -----------------------------------------------------\n"
		<< " 1. Eb/N0 [dB]						\n"
		<< " 2. Es/N0 [dB]						\n"
		<< " 3. BER								\n"
		<< " 4. BLER							\n"
		<< " 5. Avg. STE / info. bit			\n"
		<< " 6. Avg. COM / info. bit			\n"
		<< " 7. Avg. Goal Node / decoding		\n"
		<< " 5-1. Avg. Binary STE / info. bit   \n"
		<< " 8-1. Worst Case, Avg. STE / info.bit	\n"
		<< " 8-2. Worst Case, Avg. COM / info.bit	\n"
		<< " 9. Worst Case, Goal Node			\n"
		<< " 9-1. Worst_OSC_Ratio				\n"
		<< "10. Max Used STACK					\n"
		<< "--------------------------------	\n"
		<< "11. DM_STE / info. bit				\n"
		<< "12. DM_COM / info. bit				\n"
		<< "13. DM_DN / decoding				\n"
		<< "--------------------------------	\n"
		<< "14. Spending Time[s]				\n"
		<< "15. Decoding Rate[ block / s ]		\n" 
		<< "-----------------------------------------------------\n"
		<< endl;
	if (key_in_data.Eb_N0_1_or_Es_N0_2 == 1) {
		fp
		<< "Eb/N0 Test	            			\n";
	}
	else if(key_in_data.Eb_N0_1_or_Es_N0_2 == 2) {
		fp
		<< "Es/N0 Test	            			\n";
	}
	fp.close();
};

void ClearFile(string name)
{
	using namespace std;
	string filename(name);
	fstream fp;
	fp.open(filename, ios::out | ios::trunc);
	if (!fp)
		cout << "Fail to open file: " << filename << endl;
	fp.close();
};

void WriteFile(string name, double input)
{
	using namespace std;
	string filename(name);
	fstream fp;
	fp.open(filename, ios::out | ios::app);
	if (!fp) cout << "Fail to open file: " << filename << endl;
	fp << input << "\n";
	fp.close();
};


void WriteFile(string name, string input)
{
	using namespace std;
	string filename(name);
	fstream fp;
	fp.open(filename, ios::out | ios::app);
	if (!fp) cout << "Fail to open file: " << filename << endl;

	fp << input << "\n";
	fp.close();
};

void WriteFile(string name, vector<size_t> input)
{
	using namespace std;
	string filename(name);
	fstream fp;
	fp.open(filename, ios::out | ios::app);
	if (!fp) cout << "Fail to open file: " << filename << endl;
	
	for (size_t i(0); i < input.size(); i++)
		if (input[i] == 1) fp << i << " ";
	fp << "\n";
	fp.close();
}

/*
void WriteFile2(string name, vector<size_t> input)
{
	using namespace std;
	string filename(name);
	fstream fp;
	fp.open(filename, ios::out | ios::trunc);
	if (!fp) cout << "Fail to open file: " << filename << endl;

	for (size_t i(0); i < input.size(); i++)
		fp << input[i] << " ";
	fp << "\n";
	fp.close();
}*/

