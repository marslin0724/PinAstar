#include "PinAstar/LinearBlockCodes.h"
#include "stdio.h"
/*
enum OpType {
	OT_ADD,
};

template<class T1, class T2>
class VecSum {
	OpType type_;
	const T1& op1_;
	const T2& op2_;
public:
	VecSum(OpType type, const T1& op1, const T2& op2) : type_(type), op1_(op1), op2_(op2) {}
	
	unsigned int operator[](const unsigned int i) const {
		switch (type_) {
		case OT_ADD: return op1_[i] ^ op2_[i];
		default: throw "bad type";
		}
	}
};

template<class T1, class T2>
VecSum<T1, T2> operator ^ (const T1& t1, const T2& t2) {
	return VecSum<T1, T2>(OT_ADD, t1, t2);
}*/


// Linear block code 
/******************************************************************/
size_t G_QR24_12[1] = { 06165000000 };
size_t G_QR48_24[1] = { 04307335700 };
size_t G_QR80_40[2] = { 01143571726, 04067000000 };
size_t G_QR90_45[2] = { 05523367773, 06625500000 };
size_t G_QR98_49[2] = { 01744155434, 03330237000 };
size_t G_QR104_52[2] = { 01307024764, 00757141300 };
size_t G_QR192_96[4] = { 07741633224, 03422557065, 06354432042, 00100000000 };                                 // 政晏


size_t G_BCH128_64[3] = { 01206534025, 05707731000, 04500000000 };
size_t G_BCH128_78[2] = { 05446000435, 04260232000};  // Chou-Yin Code Rate=0.6093 
size_t G_BCH128_99[1] = { 03447023271 };
size_t G_BCH64_36[1] = { 04156402114 };
size_t G_BCH256_45[8] = { 01520205605, 05234161131, 01013463764, 02370156367, 00024470762, 03730332021, 05702505154, 0100000000 };
size_t G_BCH256_99[6] = { 01065666725, 03473174222, 07414162015, 07433225241, 01076432303, 04310000000 };
size_t G_BCH256_131[5] = { 02157133314, 07151015126, 01250277442, 01420241654, 07100000000 };
//Pin
size_t G_BCH128_36[4] = { 06314171555, 02441721117, 05137164367, 02000000000 };

// generator polynomial
vector<__int8> poly_4 = { 1,0,0,1,1 };

// RM vectors
vector<__int8> RM_v0 = { 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1 };
vector<__int8> RM_v1 = { 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1 };

// 張家輔學長 博論 例子
__int8 G_example[4][8] =
{
{1,0,0,0,1,0,1,1},
{0,1,0,0,1,1,1,0},
{0,0,1,0,1,1,0,1},
{0,0,0,1,0,1,1,1}
};

//(7,4) hamming code
__int8 G_hamming_7_4[4][7] =
{
{1,0,0,0,1,0,1},
{0,1,0,0,1,1,1},
{0,0,1,0,1,1,0},
{0,0,0,1,0,1,1}
};

/******************************************************************/

//#define Polar_AStar_Combination     102
/******************************************************************/

void CODE::ChooseCode()
{
	cout
		<< "\n Binary Liear Block Codes : \n"
		<< " \n**** Golay Code **** \n\n"
		<< " " << Golay_24_12 << ": (24,12,8) Golay Code \n"

		<< " \n****** QR Code ***** \n\n"
		<< " " << QR_48_24 << ": (48,24,12) QR Code	\n"
		<< " " << QR_80_40 << ": (80,40,16) QR Code	\n"
		<< " " << QR_90_45 << ": (90,45,18) QR Code	\n"
		<< " " << QR_98_49 << ": (98,49,16) QR Code	\n"
		<< " " << QR_104_52 << ": (104,52,20) QR Code	\n"

		<< " \n****** BCH Code ***** \n\n"
		<< " " << BCH_128_64 << ": (128,64) BCH Code		\n"
		<< " " << BCH_128_99 << ": (128,99) BCH Code		\n"
		<< " " << BCH_256_45 << ": (256,45) BCH Code		\n"
		<< " " << BCH_256_99 << ": (256,99) BCH Code		\n"
		<< " " << BCH_256_131 << ": (256,131) BCH Code	\n"
		<< " " << BCH_128_36 << ": (128,36) BCH Code \n"
		<< " " << BCH_64_36_Rep << ": (64,36)eBCH + repeatition\n"

		<< "\n****** RS Code ***** \n\n"
		<< " " << RSGF2_64_36 << ": (64,36) RS convert to GF(2)	\n"
		<< " " << RS_16_9_concat_8_4 << ": RS(16,9,8) concatenated (8,4)HM \n"
		<< " " << HM_8_4 << ": Hamming 8 4 code \n"
		<< " " << HM2_1  << ": Hamming 2 1 code \n"
		<< " " << RS_16_9_concat_8_4_interleaver << ": RS16,9 + (8,4) + interleaver \n"
		<< " " << RSGF2_64_36_Rep << ": RSGF2 (64,36) Repeatition\n"
		<< " " << RS16_15		<< ": RS16_15\n"
		<< " " << RS16_15_concate_HM8_4 << ": RS16_15 concatenate HM8_4\n"
		<< " " << RS16_14 << ": RS16_14\n"
		<< " " << RS16_14_concate_HM8_4 << ": RS16_14 concatenate HM8_4\n"
		<< " " << RS16_13 << ": RS16_13\n"
		<< " " << RS16_13_concate_HM8_4 << ": RS16_13 concatenate HM8_4\n"
		<< " " << RS_16_9_concat_8_4_new << ": RS_16_9_concat_8_4_new\n "
		<< " " << C6_5					 << ": (6,5,2) code\n"
		<< " " << RS_32_19_concat_6_5   << ": RS(32,19) + (6,5,2)\n"
		<< " " << RSconcate_192_96		<< ": RSconcate_192_96\n"
		<< " " << RS16_9_concat_16_8	<< ": RS16_9_concat_16_8\n"
		<< " " << RS16_9_concat_24_12	<< ": RS16_9_concat_24_12\n"
		<< " " << RS16_8_concat_16_8	<< ": RS16_8_concat_16_8\n"
		<< " " << RS16_10_concat_16_8	<< ": RS16_10_concat_16_8\n"
		<< " " << BCH_64_36				<< ": BCH_64_36\n"
		<< " " << RS16_10_concat_8_4	<< ": RS16_10_concat_8_4\n"

		<< " \n****** Others ***** \n\n"
		<< " " << LDPC_96_48 << ": (96,48) LDPC			\n"
		<< " " << Polar_120_40 << ": (120,40) Polar Code	\n"
		<< " " << ReadFromTxt << ": Read From .txt		\n"
		<< " " << RandomCode << ": Random Code			\n"
		<< " " << BPSK << ": BPSK, Uncoded			\n"
		<< " " << LDPC_16_8 << ": (16,8) code \n"

		//<< " \n**** Integration *** \n\n"
		//<< " " << BCH_128_LDPC_1024 << ": six (128,99) BCH Code with LDPC (1024,768)			\n"
		//<< " " << BCH_128_78_ProductCode << ": Product code of { (128,78,16) BCH Code & (8,7,2) parity check code }\n"
		//<< " " << BCH_128_78_ProductCode_Rev << ": Product code of { (8,7,2) parity check code & (128,78,16) BCH Code }\n"

		<< " \n******* Test ******* \n\n"
		<< " " << BCH_128_78 << ": (128,78) BCH Code		     \n"
		<< " " << QR_192_96 << ": (196,92) QR Code		         \n"
		<< " " << RM_16_11 << ": (16,11) RM Code -> (r,m) = (2,4)\n"
		<< " " << RM_32_6 << ": (32,6) RM Code -> (r,m) = (1,5)\n"
		<< " " << RM_64_42 << ": (64,42) RM Code -> (r,m) = (3,6)\n"
		<< " " << Hamming_Code_15_11 << ": (15,11) Hamming Code  \n"
		<< " \n******* PoHan Tests ******* \n\n"
		<< " " << double_RandomCode << ": double Random Code			\n"
		<< " " << Parity_Check << ": Parity Check Code			\n"
		<< " " << Parity_CheckII << ": Parity Check Code Hybrid			\n"
		<< " " << Parity_CheckIII << ": Parity Check CodeIII			\n"
		<< " " << Random_ParityCheck << ": Random Position Parity Check Code			\n"
		<< " " << ReedSolomonConcatenated << ": Reed-Solomon Concatenated Code (192,96) \n"

		<< " \n******* 5G codes ******* \n\n"
		<< " " << LDPC_256_192 << ": LDPC Code (256,192)        \n"
		<< " " << LDPC_128_64  << ": LDPC Code (128,64)         \n"
		<< " " << LDPC_192_96  << ": LDPC Code (192,96)         \n"
		<< " " << LDPC_256_128 << ": LDPC Code (256,128)        \n"
		<< " " << LDPC_512_256 << ": LDPC Code (512,256)        \n"
		<< " " << LDPC_224_128 << ": LDPC Code (224,128)        \n"
		<< " " << LDPC_160_128 << ": LDPC Code (160,128)        \n"
		<< " " << LDPC_192_128 << ": LDPC Code (192,128)        \n"
		<< " " << Polar_Code   << ": Polar Code                 \n"

		<< " \n**** Integration & Some complicated codes *** \n\n"
		<< " " << Polar_256_192_AStar_192_128   << ": Polar Code (256,192) contenate RC (192,128) \n"
		<< " " << LDPC_256_192_AStar_192_128    << ": LDPC Code (256,192) contenate RC (192,128) \n"
		<< " " << Hamming_255_187_Astar_187_128 << ": Outer RC (187,128) Inner 17 Hamming Code (255,187)\n"
		<< " " << RM_256_192_Astar_192_128      << ": Outer RC (176,128) Inner 16 RM Code (16,11), Total:(256,176)\n"
		<< " " << RM_256_168_Astar_168_128      << ": Outer RC (168,128) Inner  4 RM Code (64,42), Total:(256,168)\n"
		<< " " << Turbo_256_128_with_2LDPC      << ": (256,128) Turbo Code with 2 (192,128) LDPC Code\n"
		<< " " << Turbo_256_128_with_2LDPC_Punc << ": (256,128) Turbo Code with 2 (192,128) LDPC Code, Punctured version\n"
		<< " " << Turbo_256_128_with_2LDPC_ver2 << ": (256,128) Turbo Code with (160,128) LDPC Code & (224,128) LDPC Code\n"
		<< " " << Turbo_320_128_with_2LDPC      << ": (320,128) Turbo Code with (224,128) LDPC Code & (224,128) LDPC Code\n"
		<< " " << RandomCode_with_LDPC          << ": RandomCode using SPA\n"
		<< " " << Turbo_256_128_RC_LDPC_192_128 << ": (256,128) Turbo Code with RC & LDPC Code(192, 128)\n"
		<< " " << Polar_256_128_Astar           << ": (256,128) Polar Code with Astar \n"
		<< " " << Polar_with_Astar_middle_index << ": (256,128 + m) Polar Code with Outer Code (2m, m) \n"
		
		<< "\n Select a number : ";
	Code_number = 0;
	std::cin >> Code_number;
	Inner_Col_Number = 0;
	Short_Inner_Col_Number = 0;
	Punctured_Number = 0;
	LDPC_FUNC LDPC;
	char str[10];
	size_t paritybitsperrow(0), paritybitsperrowfirst(0), paritybitsperrowsecond(0);
	//interleaver
	vector<int> tmp{ 0,17,36,53,1,20,37,56,4,21,40,57,5,24,41,58,8,25,44,59,9,28,45,60,12,
		29,48,61,13,32,49,62,16,33,52,63,2,15,30,43
		,3,18,31,46,6,19,34,47,7,22,35,50,10,23,38,51,11,26,39,54,14,27,42,55 };
	switch (Code_number)
	{
	case Golay_24_12://Golay code   24, 12, 8
		Col_Number = 24;
		Row_Number = 12;
		Title = "(24,12,8) Golay Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_QR24_12, G);
		break;

	case QR_48_24://QR code  48, 24, 12
		Col_Number = 48;
		Row_Number = 24;
		Title = "(48,24,12) QR Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_QR48_24, G);
		break;

	case QR_80_40://QR code   80, 40, 16
		Col_Number = 80;
		Row_Number = 40;
		Title = "(80,40,16) QR Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_QR80_40, G);
		break;

	case QR_90_45://QR code    90, 45, 18
		Col_Number = 90;
		Row_Number = 45;
		Title = "(90,45,18) QR Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_QR90_45, G);
		break;

	case QR_98_49://QR code   98, 49, 16
		Col_Number = 98;
		Row_Number = 49;
		Title = "(98,49,16) QR Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_QR98_49, G);
		break;

	case QR_104_52://QR code   104, 52, 20
		Col_Number = 104;
		Row_Number = 52;
		Title = "(104,52,20) QR Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_QR104_52, G);
		break;

	case Polar_120_40://Polar code 120,40
		Title = "(120,40) Polar Code";
		ReadFile_GeneratorMatrix("Polar_120_40.txt",G);
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		MatrixForm_to_Generator(G);
		break;

	case LDPC_96_48: //LDPC 96,48
		Col_Number = 96;
		Row_Number = 48;
		Title = "(96,48) LDPC ";
		ReadFile_GeneratorMatrix("LDPC_96_48.txt", G);
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		MatrixForm_to_Generator(G);
		break;

	case BCH_128_64://BCH code  128, 64, 22
		Col_Number = 128;
		Row_Number = 64;
		Title = "(128,64,22) BCH Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH128_64, G);
		Convert_SystematicG_to_H(H, G);
		G._H = H._matrix;
		break;

	case BCH_128_99://BCH code  128, 99, 22
		Col_Number = 128;
		Row_Number = 99;
		Title = "(128,99) BCH Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH128_99, G);
		break;

	case BCH_256_45://BCH code  
		Col_Number = 256;
		Row_Number = 45;
		Title = "(256,45) BCH Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH256_45, G);
		break;

	case BCH_256_99:
		Col_Number = 256;
		Row_Number = 99;
		Title = "(256,99) BCH Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH256_99, G);
		break;

	case BCH_256_131: //BCH code 256,131
		Col_Number = 256;
		Row_Number = 131;
		Title = "(256,131) BCH Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH256_131, G);
		break;
	case BCH_128_36: //BCH code 256,131
		Col_Number = 128;
		Row_Number = 36;
		Title = "(128,36) BCH Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH128_36, G);
		Convert_SystematicG_to_H(H, G);
		G._H = H._matrix;
		break;
	case BCH_64_36_Rep:
		Col_Number = 64;
		Row_Number = 36;
		Title = "(64,36) BCH Code + repeatition";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH64_36, G);
		G._matrix_inner = G._matrix;
		repeatition(G);
		G._matrix_outer.resize(1);
		G._matrix_outer.at(0).resize(2);
		G._matrix_outer.at(0).at(0) = G._matrix_outer.at(0).at(1) = 1;
		Col_Number = 128;
		break;
	case BCH_64_36://BCH code  
		Col_Number = 64;
		Row_Number = 36;
		Title = "(64,36) BCH Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH64_36, G);
		break;
	case RSGF2_64_36:
		Title = "(64,36)Reed-Solomon Code ";
		ReadFile_GeneratorMatrix("RSGF2_64_36.txt", G);
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		MatrixForm_to_Generator(G);
		Convert_SystematicG_to_H(H, G);
		G._H = H._matrix;
		G._matrix_outer = H._matrix;
		break;
	case RS_16_9_concat_8_4:
		Title = "RS_16_9_concat_8_4";
		ReadFile_GeneratorMatrix("RSGF2_64_36.txt", G_Inner);
		ReadFile_GeneratorMatrix("HM8_4_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		H.Building_Empty_Matrix(4, 8);
		ReadFile_GeneratorMatrix("HM8_4.txt", H);//讀取H矩陣
		G._matrix_outer = H._matrix;
		//RS16_9 parity check
		Convert_SystematicG_to_H(H, G_Inner);
		G._H = H._matrix;
		ReadFile_GeneratorMatrix("G8_4.txt", G_);
		G.subG = G_._matrix;
		//test
		//GH_test(G_Inner, H);
		//
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case HM_8_4:
		Title = "Hamming_8_4";
		ReadFile_GeneratorMatrix("G8_4.txt",G);
		H.Building_Empty_Matrix(4, 8);
		ReadFile_GeneratorMatrix("HM8_4.txt", H);//讀取H矩陣
		G._matrix_outer = H._matrix;
		G._H = H._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case HM2_1:
		Title = "Hamming_2_1";
		G.Building_Empty_Matrix(1, 2);
		G._matrix.at(0).at(0) = G._matrix.at(0).at(1) = 1;
		G._matrix_outer = G._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS_16_9_concat_8_4_interleaver:
		Title = "RS_16_9_concat_8_4_interleaver";
		ReadFile_GeneratorMatrix("RSGF2_64_36.txt", G_Inner);
		ReadFile_GeneratorMatrix("HM8_4_concate_interleaver.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		H.Building_Empty_Matrix(4, 8);
		ReadFile_GeneratorMatrix("HM8_4.txt", H);//讀取H矩陣 
		G.interleaver = tmp;
		G._matrix_outer = H._matrix;
		ReadFile_GeneratorMatrix("G8_4.txt", G_);
		G.subG = G_._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RSGF2_64_36_Rep:
		Title = "(64,36)Reed-Solomon Code ";
		ReadFile_GeneratorMatrix("RSGF2_64_36.txt", G);
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		MatrixForm_to_Generator(G);
		repeatition(G);
		Col_Number *= 2;
		break;
	case LDPC_16_8:
		Title = "LDPC_16_8";
		ReadFile_GeneratorMatrix("LDPC16_8.txt", G);
		H.Building_Empty_Matrix(8, 16);
		Convert_SystematicG_to_H(H, G);
		G._H = H._matrix;
		G._matrix_outer = H._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS16_15:
		Title = "(16,15)Reed-Solomon Code ";
		ReadFile_GeneratorMatrix("RS16_15.txt", G);
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		MatrixForm_to_Generator(G);
		Convert_SystematicG_to_H(H, G);
		G._H = H._matrix;
		G._matrix_outer = H._matrix;
		break;
	case RS16_15_concate_HM8_4:
		Title = "RS_16_15_concat_8_4";
		ReadFile_GeneratorMatrix("RS16_15.txt", G_Inner);
		ReadFile_GeneratorMatrix("HM8_4_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		H.Building_Empty_Matrix(4, 8);
		ReadFile_GeneratorMatrix("HM8_4.txt", H);//讀取H矩陣
		G._matrix_outer = H._matrix;
		ReadFile_GeneratorMatrix("G8_4.txt", G_);
		G.subG = G_._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS16_14:
		Title = "(16,14)Reed-Solomon Code ";
		ReadFile_GeneratorMatrix("RS16_14.txt", G);
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		MatrixForm_to_Generator(G);
		Convert_SystematicG_to_H(H, G);
		ReadFile_GeneratorMatrix("HM8_4.txt", G_);
		G.subG = G_._matrix;
		G._H = H._matrix;
		G._matrix_outer = H._matrix;
		break;
	case RS16_14_concate_HM8_4:
		Title = "RS_16_14_concat_8_4";
		ReadFile_GeneratorMatrix("RS16_14.txt", G_Inner);
		ReadFile_GeneratorMatrix("HM8_4_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		H.Building_Empty_Matrix(4, 8);
		ReadFile_GeneratorMatrix("HM8_4.txt", H);//讀取H矩陣
		G._matrix_outer = H._matrix;
		ReadFile_GeneratorMatrix("G8_4.txt", G_);
		G.subG = G_._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS16_13:
		Title = "(16,13)Reed-Solomon Code ";
		ReadFile_GeneratorMatrix("RS16_13.txt", G);
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		MatrixForm_to_Generator(G);
		Convert_SystematicG_to_H(H, G);
		G._H = H._matrix;
		G._matrix_outer = H._matrix;
		break;
	case RS16_13_concate_HM8_4:
		Title = "RS_16_13_concat_8_4";
		ReadFile_GeneratorMatrix("RS16_13.txt", G_Inner);
		ReadFile_GeneratorMatrix("HM8_4_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		H.Building_Empty_Matrix(4, 8);
		ReadFile_GeneratorMatrix("HM8_4.txt", H);//讀取H矩陣
		G._matrix_outer = H._matrix;
		ReadFile_GeneratorMatrix("G8_4.txt", G_);
		G.subG = G_._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS_16_9_concat_8_4_new:
		Title = "RS_16_9_concat_8_4_new";
		H.Building_Empty_Matrix(28, 64);
		ReadFile_GeneratorMatrix("RS16_9_H.txt", H);
		MatrixForm_to_Generator(H);
		exchang_column(H, 28);
		G_Inner.Building_Empty_Matrix(36, 64);
		convert_H_G(H, G_Inner);
		//test
		GH_test(G_Inner, H);
		//
		ReadFile_GeneratorMatrix("HM8_4_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		G._H = H._matrix;
		H.Building_Empty_Matrix(4, 8);
		ReadFile_GeneratorMatrix("HM8_4.txt", H);//讀取H矩陣
		G._matrix_outer = H._matrix;
		ReadFile_GeneratorMatrix("G8_4.txt", G_);
		G.subG = G_._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case C6_5:
		Title = "(6,5) code";
		ReadFile_GeneratorMatrix("G6_5.txt", G);
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS_32_19_concat_6_5:
		Title = "RS_32_19_concat_6_5";
		ReadFile_GeneratorMatrix("RS32_19.txt", G_Inner);
		ReadFile_GeneratorMatrix("C6_5_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		H.Building_Empty_Matrix(4, 8);
		//RS32_19 parity check
		Convert_SystematicG_to_H(H, G_Inner);
		G._H = H._matrix;
		ReadFile_GeneratorMatrix("G6_5.txt", G_);
		G.subG = G_._matrix;
		//test
		//GH_test(G_Inner, H);
		//
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RSconcate_192_96:
		Title = "RSconcate_192_96";
		ReadFile_GeneratorMatrix("RS16_9.txt", G_Inner);
		ReadFile_GeneratorMatrix("C6_5_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G.Row_number++;
		G._matrix.resize(G.Col_number);
		G._matrix.at(G.Row_number - 1) = vector<__int8>(G.Col_number, 1);
		MatrixForm_to_Generator(G);
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS16_9_concat_16_8:
		Title = "RS16_9_concat_16_8";
		ReadFile_GeneratorMatrix("RSGF2_64_36.txt", G_Inner);
		ReadFile_GeneratorMatrix("C16_8_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		H.Building_Empty_Matrix(4, 8);
		//RS32_19 parity check
		Convert_SystematicG_to_H(H, G_Inner);
		G._H = H._matrix;
		ReadFile_GeneratorMatrix("LDPC16_8.txt", G_);
		G.subG = G_._matrix;
		//test
		//GH_test(G_Inner, H);
		//
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS16_9_concat_24_12:
		Title = "RS16_9_concat_24_12";
		ReadFile_GeneratorMatrix("RSGF2_64_36.txt", G_Inner);
		ReadFile_GeneratorMatrix("C24_12_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		H.Building_Empty_Matrix(4, 8);
		//RS32_19 parity check
		Convert_SystematicG_to_H(H, G_Inner);
		G._H = H._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS16_8_concat_16_8:
		Title = "RS16_8_concat_16_8";
		ReadFile_GeneratorMatrix("RS16_8.txt", G_Inner);
		ReadFile_GeneratorMatrix("C16_8_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		Convert_SystematicG_to_H(H, G_Inner);
		G._H = H._matrix;
		ReadFile_GeneratorMatrix("LDPC16_8.txt", G_);
		G.subG = G_._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS16_10_concat_16_8:
		Title = "RS16_10_concat_16_8";
		ReadFile_GeneratorMatrix("RS16_10.txt", G_Inner);
		ReadFile_GeneratorMatrix("C16_8_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		Convert_SystematicG_to_H(H, G_Inner);
		G._H = H._matrix;
		ReadFile_GeneratorMatrix("LDPC16_8.txt", G_);
		G.subG = G_._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS16_11_concat_16_8:
		Title = "RS16_11_concat_16_8";
		ReadFile_GeneratorMatrix("RS16_11.txt", G_Inner);
		ReadFile_GeneratorMatrix("C16_8_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		Convert_SystematicG_to_H(H, G_Inner);
		G._H = H._matrix;
		ReadFile_GeneratorMatrix("LDPC16_8.txt", G_);
		G.subG = G_._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case RS16_10_concat_8_4:
		Title = "RS16_10_concat_8_4";
		ReadFile_GeneratorMatrix("RS16_10.txt", G_Inner);
		ReadFile_GeneratorMatrix("HM8_4_concatenate.txt", G_);
		G.Building_Empty_Matrix(G_Inner.Row_number, G_.Col_number);
		Matrix_Mul(G_Inner, G_, G);
		G._matrix_inner = G_Inner._matrix;
		Convert_SystematicG_to_H(H, G_Inner);
		G._H = H._matrix;
		ReadFile_GeneratorMatrix("G8_4.txt", G_);
		G.subG = G_._matrix;
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		break;
	case ReadFromTxt:
		Title = "Read From .txt";
		ReadFile_GeneratorMatrix("Other_G.txt", G);
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		MatrixForm_to_Generator(G);

		sprintf_s(str, sizeof(str),"%d", Col_Number);
		Title = Title + ",(" + str + ",";
		sprintf_s(str,sizeof(str), "%d", Row_Number);
		Title = Title + str + ") ";
		break;
	case RandomCode:
		cout << "\n Codeword Length : ";
		cin >> Col_Number;
		cout << "\n Message Length: ";
		cin >> Row_Number;
		Title = "Random Code ( P(bit=1)=P(bit=0)=0.5 )";
		//Yin
		if (CRC_Check == TRUE) {
			Row_Number += CRC_length;
		}
		//Yin
		sprintf_s(str, sizeof(str),"%d", Col_Number);
		Title = Title + ",(" + str + ",";
		sprintf_s(str, sizeof(str), "%d", Row_Number);
		Title = Title + str + ") ";
		CreateRandomCode(Row_Number, Col_Number);
		ReadFile_GeneratorMatrix("RandomCode.txt", G);
		MatrixForm_to_Generator(G);
		break;
	
	case BPSK:
		Row_Number = 1;
		Col_Number = 1;
		Title = "BPSK , Uncoded Signalling";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		G._matrix[0][0] = 1;
		break;
	/*case BCH_128_LDPC_1024: // BCH code: 128,99,22
		Col_Number = 128;
		Row_Number = 99;
		Title = "six (128,99) BCH Code with LDPC (1024,768)";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH128_99, G);
		break;
	case BCH_128_78_ProductCode:
		Col_Number = 128;
		Row_Number = 78;
		Title = "Product code of { (128,78,16) BCH Code & (8,7,2) parity check code }";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH128_78, G);
		break;
	case BCH_128_78_ProductCode_Rev:
		Col_Number = 128;
		Row_Number = 78;
		Title = "Product code of { (8,7,2) parity check code & (128,78,16) BCH Code }";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH128_78, G);
		break;*/

	case BCH_128_78://BCH code  128, 78
		Col_Number = 128;
		Row_Number = 78;
		Title = "(128,78) BCH Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_BCH128_78, G);
		break;
	case QR_192_96://QR code  192, 96
		Col_Number = 192;
		Row_Number = 96;
		Title = "(192,96) QR Code";
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		Extended_BCH_Generator_Matrix(G_QR192_96, G);
		break;
	case RM_16_11:
		Col_Number = 16;   // 16 
		Row_Number = 11;   // 11 
		Title = "(16,11) RM code a.k.a. (r,m) = (2,4)";
		RM_16_11_Generator_Matrix(G);
		break;
	case RM_32_6:
		Col_Number = 32;   // 32 
		Row_Number = 6;    // 6 
		Title = "(32,6) RM code a.k.a. (r,m) = (1,5)";
		RM_32_6_Generator_Matrix(G);
		break;
	case RM_64_42:
		Col_Number = 64;   // 64 
		Row_Number = 42;   // 42 
		Title = "(64,42) RM code a.k.a. (r,m) = (3,6)";
		RM_order3_Generator_Matrix(6, G, G_);
		break;
	case Hamming_Code_15_11:
		Col_Number = 15;   // 15 
		Row_Number = 11;   // 11 
		Title = "(15,11) Hamming code";
		HammingCode_Generator_Matrix(G, H, poly_4, 15);
		break;
	case LDPC_16_12:
		Col_Number = 16;   // 16 
		Row_Number = 12;   // 12 
		Title = "(16,12) LDPC code";
		H.Building_Empty_Matrix(Col_Number - Row_Number, Col_Number);
		LDPC_Code_16_12(H);
		Permutation_Seq.resize(Col_Number, 0);
		LDPC.H_G_convertor_G_P_I(H, Permutation_Seq, G);
		break;
	//PoHan
	case double_RandomCode:
	{
		size_t different_probability_row_position;
		double probability1, probability2;
		cout << "\n Codeword Length: ";
		cin >> Col_Number;
		cout << "\n Message Length: ";
		cin >> Row_Number;
		cout << "\n Different Probability Row: ";
		cin >> different_probability_row_position;
		cout << "\n Probability 1: ";
		cin >> probability1;
		cout << "\n Probability 2: ";
		cin >> probability2;
		H.Col_number = Col_Number;
		H.Row_number = Col_Number - Row_Number;
		G.Col_number = Col_Number;
		G.Row_number = Row_Number;
		H.Building_Empty_Matrix(H.Row_number, H.Col_number);
		CreateRandomCode(0, different_probability_row_position, 0.1, H);
		CreateRandomCode(different_probability_row_position, Row_Number, 0.75, H);
		Permutation_Seq.resize(Col_Number, 0);
		H_G_convertor_G_P_I(H, Permutation_Seq, G);
	}
		break;
	case Parity_Check:
		cout << "\n Codeword Length: ";
		cin >> Col_Number;
		cout << "\n Message Length: ";
		cin >> Row_Number;
		cout << "\n Parity Bits Per Row: ";
		cin >> paritybitsperrow;
		H.Col_number = Col_Number;
		H.Row_number = Col_Number - Row_Number;
		G.Col_number = Col_Number;
		G.Row_number = Row_Number;
		H.Building_Empty_Matrix(H.Row_number, H.Col_number);

		for (size_t i(0); i < H.Col_number / paritybitsperrow; i++) {
			for (size_t j(0); j < paritybitsperrow; j++) {
				H._matrix[i][paritybitsperrow * i + j] = 1;
			}
		}
			
		FillRandomBits(H.Col_number / paritybitsperrow, H);
		Permutation_Seq.resize(Col_Number, 0);
		H_G_convertor_G_P_I(H, Permutation_Seq, G);
		H.Low_Density_Row_Index.resize(H.Col_number / paritybitsperrow, 0);
		for (size_t i(0); i < H.Row_number; i++) {
			size_t k(0);
			for (size_t j(0); j < H.Col_number; j++) {
				if ((int)H._matrix[i][j] == 1) {
					k++;
				}
			}
			if (k == paritybitsperrow) {
				H.Low_Density_Row_Index.at(H.Num_Low_Density_Row) = i;
				H.Num_Low_Density_Row++;
			}
		}
	break;

	case Parity_CheckII: // 每列bits數 一半4個一半8個
	{
		cout << "\n Codeword Length: ";
		cin >> Col_Number;
		cout << "\n Message Length: ";
		cin >> Row_Number;
		cout << "\n Parity Bits Per Row for first half : ";
		cin >> paritybitsperrowfirst;
		cout << "\n Parity Bits Per Row for second half : ";
		cin >> paritybitsperrowsecond;
		H.Col_number = Col_Number;
		H.Row_number = Col_Number - Row_Number;
		G.Col_number = Col_Number;
		G.Row_number = Row_Number;
		H.Building_Empty_Matrix(H.Row_number, H.Col_number);
		
		for (size_t i(0); i < H.Col_number / (paritybitsperrowfirst *2); i++) {
			H._matrix[i][paritybitsperrowfirst * i] = 1;
			H._matrix[i][paritybitsperrowfirst * i + 1] = 1;
			H._matrix[i][paritybitsperrowfirst * i + 2] = 1;
			H._matrix[i][paritybitsperrowfirst * i + 3] = 1;
		}
		for (size_t i(0); i < H.Col_number / (paritybitsperrowsecond * 2); i++) {
			H._matrix[i + H.Col_number / (paritybitsperrowfirst * 2)][paritybitsperrowsecond * (i + H.Col_number / (paritybitsperrowfirst * 4))] = 1;
			H._matrix[i + H.Col_number / (paritybitsperrowfirst * 2)][paritybitsperrowsecond * (i + H.Col_number / (paritybitsperrowfirst * 4)) + 1] = 1;
			H._matrix[i + H.Col_number / (paritybitsperrowfirst * 2)][paritybitsperrowsecond * (i + H.Col_number / (paritybitsperrowfirst * 4)) + 2] = 1;
			H._matrix[i + H.Col_number / (paritybitsperrowfirst * 2)][paritybitsperrowsecond * (i + H.Col_number / (paritybitsperrowfirst * 4)) + 3] = 1;
			H._matrix[i + H.Col_number / (paritybitsperrowfirst * 2)][paritybitsperrowsecond * (i + H.Col_number / (paritybitsperrowfirst * 4)) + 4] = 1;
			H._matrix[i + H.Col_number / (paritybitsperrowfirst * 2)][paritybitsperrowsecond * (i + H.Col_number / (paritybitsperrowfirst * 4)) + 5] = 1;
			H._matrix[i + H.Col_number / (paritybitsperrowfirst * 2)][paritybitsperrowsecond * (i + H.Col_number / (paritybitsperrowfirst * 4)) + 6] = 1;
			H._matrix[i + H.Col_number / (paritybitsperrowfirst * 2)][paritybitsperrowsecond * (i + H.Col_number / (paritybitsperrowfirst * 4)) + 7] = 1;
		}
		FillRandomBits((H.Col_number / (paritybitsperrowfirst * 2) + H.Col_number / (paritybitsperrowsecond * 2)), H);
		/*
		for (size_t i(0); i < H.Row_number; ++i) {
			for (size_t j(0); j < H.Col_number; ++j) {
				cout << (int)H._matrix[i][j];
			}
			cout << "\n";
		}
		system("PAUSE");
		*/
		Permutation_Seq.resize(Col_Number, 0);
		H_G_convertor_G_P_I(H, Permutation_Seq, G);
		H.Low_Density_Row_Index.resize(H.Col_number / (paritybitsperrowfirst * 2), 0);
		for (size_t i(0); i < H.Row_number; i++) {
			size_t k(0);
			for (size_t j(0); j < H.Col_number; j++) {
				if ((int)H._matrix[i][j] == 1) {
					k++;
				}
			}
			if (k == paritybitsperrowfirst) {
				H.Low_Density_Row_Index.at(H.Num_Low_Density_Row) = i;
				H.Num_Low_Density_Row++;
			}
		}
	}
	break;

	case Parity_CheckIII: // H interleave
	{
		cout << "\n Codeword Length: ";
		cin >> Col_Number;
		cout << "\n Message Length: ";
		cin >> Row_Number;
		cout << "\n Parity Bits Per Row: ";
		cin >> paritybitsperrow;
		H.Col_number = Col_Number;
		H.Row_number = Col_Number - Row_Number;
		G.Col_number = Col_Number;
		G.Row_number = Row_Number;
		H.Building_Empty_Matrix(H.Row_number, H.Col_number);
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		srand((unsigned int)time(NULL));
		vector <int> interleaver(H.Col_number, 0);
		int j;
		for (size_t i(1); i <= H.Col_number; i++) {
			do {
				interleaver[i - 1] = (rand() % H.Col_number);
				for (j = 1; j < i; j++) {
					if (interleaver[i - 1] == interleaver[j - 1]) {
						break;
					}
				}
			} while (j != i);
		}
		
		for (size_t i(0); i < H.Col_number / paritybitsperrow; i++) {
			for (size_t j(0); j < paritybitsperrow; j++) {
				H._matrix[i][interleaver[paritybitsperrow * i + j]] = 1;
			}
		}
		
		FillRandomBits(H.Col_number / paritybitsperrow, H);
		Permutation_Seq.resize(Col_Number, 0);
		H_G_convertor_G_P_I(H, Permutation_Seq, G);
		H.Low_Density_Row_Index.resize(H.Col_number / paritybitsperrow, 0);
		for (size_t i(0); i < H.Row_number; i++) {
			size_t k(0);
			for (size_t j(0); j < H.Col_number; j++) {
				if ((int)H._matrix[i][j] == 1) {
					k++;
				}
			}
			if (k == paritybitsperrow) {
				H.Low_Density_Row_Index.at(H.Num_Low_Density_Row) = i;
				H.Num_Low_Density_Row++;
			}
		}
		
		/*

		FillRandomBits(H.Col_number / paritybitsperrow, H);
		Permutation_Seq.resize(Col_Number, 0);
		LDPC.H_G_convertor_G_P_I(H, Permutation_Seq, G);
		MatrixForm_to_Generator(G);
		
		for (size_t i(0); i < H.Col_number / paritybitsperrow; i++) {
			H._matrix[i][interleaver[paritybitsperrow * i]] = 1;
			H._matrix[i][interleaver[paritybitsperrow * i + 1]] = 1;
			H._matrix[i][interleaver[paritybitsperrow * i + 2]] = 1;
			H._matrix[i][interleaver[paritybitsperrow * i + 3]] = 1;
		}
		FillRandomBits(H.Col_number / paritybitsperrow, H);
		
		for (int i = 0; i < H.Row_number; ++i) {
			for (int j = 0; j < H.Col_number / 2; ++j) {
				cout << (int)H._matrix[i][j];
			}
			cout << endl;
		}
		cout << endl;
		system("PAUSE");
		Permutation_Seq.resize(Col_Number, 0);
		LDPC.H_G_convertor_G_P_I(H, Permutation_Seq, G);
		MatrixForm_to_Generator(G);
		*/
	}
	break;

	case Random_ParityCheck: 
	{
		cout << "\n Codeword Length: ";
		cin >> Col_Number;
		cout << "\n Message Length: ";
		cin >> Row_Number;
		cout << "\n Parity Bits Per Row: ";
		cin >> paritybitsperrow;
		H.Col_number = Col_Number;
		H.Row_number = Col_Number - Row_Number;
		G.Col_number = Col_Number;
		G.Row_number = Row_Number;
		H.Building_Empty_Matrix(H.Row_number, H.Col_number);
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		srand((unsigned int)time(NULL));
		vector <int> interleaver(paritybitsperrow, 0);
		int j;
		
		for (size_t i(0); i < H.Col_number / paritybitsperrow; i++) {
			for (size_t k(1); k <= paritybitsperrow; k++) {
				do {
					interleaver[k - 1] = (rand() % H.Col_number);
					for (j = 1; j < k; j++) {
						if (interleaver[k - 1] == interleaver[j - 1]) {
							break;
						}
					}
				} while (j != k);
			}
			for (j = 0; j < paritybitsperrow; j++) {
				H._matrix[i][interleaver.at(j)] = 1;
			}
		}
		FillRandomBits(H.Col_number / paritybitsperrow, H);
		Permutation_Seq.resize(Col_Number, 0);
		
		H_G_convertor_G_P_I(H, Permutation_Seq, G);
		H.Low_Density_Row_Index.resize(H.Col_number / paritybitsperrow, 0);
		for (size_t i(0); i < H.Row_number; i++) {
			size_t k(0);
			for (size_t j(0); j < H.Col_number; j++) {
				if ((int)H._matrix[i][j] == 1) {
					k++;
				}
			}
			if (k == paritybitsperrow) {
				H.Low_Density_Row_Index.at(H.Num_Low_Density_Row) = i;
				H.Num_Low_Density_Row++;
			}
		}
		
		break;
	}
	case ReedSolomonConcatenated:
		Title = "Reed-Solomon Concatenated Code ";
		ReadFile_GeneratorMatrix("ReedSolomonConcatenatedCode.txt", G);
		Col_Number = G.Col_number;
		Row_Number = G.Row_number;
		MatrixForm_to_Generator(G);
		sprintf_s(str, sizeof(str), "%d", Col_Number);
		Title = Title + ",(" + str + ",";
		sprintf_s(str, sizeof(str), "%d", Row_Number);
		Title = Title + str + ") ";
		break;
		// Concatenation or Complicated Test在90以後
	case LDPC_256_192:  // (320,192) puncturing-> (256,192)
		Col_Number = 320;
		Row_Number = 192;
		G.Col_number = 320;
		G.Row_number = 192;
		Title = "(256,192) LDPC Code";
		ReadFile_H_Matrix_ver2("H_256_192.txt", H);
		Permutation_Seq.resize(Col_Number, 0);
		LDPC.H_G_convertor_G_P_I(H, Permutation_Seq, G);
		Punctured_Number = 64;
		for (int i = 0; i < G.Row_number; ++i) {    // Delete Col by Puncturing
			G._matrix.at(i).erase(G._matrix.at(i).begin(), G._matrix.at(i).begin() + (Punctured_Number));
		}
		G.Col_number -= Punctured_Number;
		Col_Number = 256;
		/*
		for (int i = 0; i < G.Row_number; ++i) {
			for (int j = 0; j < G.Col_number / 2; ++j) {
				cout << (int)G._matrix[i][j];
			}
			cout << endl;
		}
		cout << endl;
		for (int i = 0; i < G.Row_number; ++i) {
			for (int j = 0; j < G.Col_number / 2; ++j) {
				cout << (int)G._matrix[i][j + (G.Col_number / 2)];
			}
			cout << endl;
		}

		system("pause"); */
		break;

	case LDPC_128_64:  // (128,64) Code
		Col_Number = 128;
		Row_Number = 64;
		G.Col_number = 128;
		G.Row_number = 64;
		Title = "(128,64) LDPC Code";
		ReadFile_H_Matrix_ver2("H_128_64.txt", H);
		Permutation_Seq.resize(Col_Number, 0);
		LDPC.H_G_convertor_G_P_I(H, Permutation_Seq, G);
		
		/*
		for (int i = 0; i < G.Row_number; ++i) {
			for (int j = 0; j < G.Col_number / 2; ++j) {
				cout << (int)G._matrix[i][j];
			}
			cout << endl;
		}
		cout << endl;
		for (int i = 0; i < G.Row_number; ++i) {
			for (int j = 0; j < G.Col_number / 2; ++j) {
				cout << (int)G._matrix[i][j + (G.Col_number / 2)];
			}
			cout << endl;
		}

		system("pause"); */
		break;
	case LDPC_192_96: //LDPC 192 96
		Col_Number = 192;
		Row_Number = 96;
		Title = "(192,96) LDPC ";
		ReadFile_GeneratorMatrix("H_192_96_z12_20200219_42.txt", H);
		Permutation_Seq.resize(Col_Number, 0);
		LDPC.H_G_convertor_G_P_I(H, Permutation_Seq, G);
		break;
	case LDPC_256_128:  // (256,128) Code
		Col_Number = 256;
		Row_Number = 128;
		G.Col_number = 256;
		G.Row_number = 128;
		Title = "(256,128) LDPC Code";
		ReadFile_H_Matrix_ver2("H_256_128.txt", H);
		Permutation_Seq.resize(Col_Number, 0);
		LDPC.H_G_convertor_G_P_I(H, Permutation_Seq, G);
		break;
	case LDPC_512_256:  // (256,128) Code
		Col_Number = 512;
		Row_Number = 256;
		G.Col_number = 512;
		G.Row_number = 256;
		Title = "(512,256) LDPC Code";
		ReadFile_H_Matrix_ver2("H_512_256.txt", H);
		Permutation_Seq.resize(Col_Number, 0);
		LDPC.H_G_convertor_G_P_I(H, Permutation_Seq, G);
		break;
	case LDPC_224_128:
		Col_Number = 224;
		Row_Number = 128;
		G.Col_number = 224;
		G.Row_number = 128;
		Title = "(224,128) LDPC Code";
		ReadFile_GeneratorMatrix("QC_LDPC_H_224_96_hm_7_3.txt", H);
		Permutation_Seq.resize(Col_Number, 0);
		LDPC.fullrank = TRUE;
		LDPC.H_G_convertor_G_I_P_ver2(H, Permutation_Seq, G);
		break;
	case LDPC_160_128:
		Col_Number = 160;
		Row_Number = 128;
		G.Col_number = 160;
		G.Row_number = 128;
		Title = "(160,128) LDPC Code";
		ReadFile_GeneratorMatrix("QC_LDPC_H_160_32_hm_10_2.txt", H);
		Permutation_Seq.resize(Col_Number, 0);
		LDPC.fullrank = FALSE;
		LDPC.H_G_convertor_G_I_P_ver2(H, Permutation_Seq, G);
		break;
	case LDPC_192_128:
		Col_Number = 192;
		Row_Number = 128;
		G.Col_number = 192;
		G.Row_number = 128;
		Title = "(192,128) LDPC Code";
		// QC_LDPC_H_192_64_hm_12_4.txt
		// QC_LDPC_H_189_63_hm_9_3.txt
		// QC_LDPC_192_128.txt
		
		LDPC.fullrank = 1;
		ReadFile_GeneratorMatrix("QC_LDPC_192_128.txt", H);
		Permutation_Seq.resize(Col_Number, 0);
		//LDPC.fullrank = FALSE;
		LDPC.H_G_convertor_G_I_P_ver2(H, Permutation_Seq, G);
		break;

	case Polar_Code:
		cout << "\n Message Length: ";
		cin >> Row_Number;
		cout << "\n Codeword Length : ";
		cin >> Col_Number;
		Title = "Polar Code ";
		sprintf_s(str, sizeof(str), "%d", Col_Number);
		Title = Title + ",(" + str + ",";
		sprintf_s(str, sizeof(str), "%d", Row_Number);
		Title = Title + str + ") ";
		G_.Building_Empty_Matrix(Col_Number, Row_Number);
		Polar_Code_G_Generator(G_, frozen, non_frozen, Channel_Parameter);
		G.Building_Empty_Matrix(Row_Number, Col_Number);
		for (int i = 0; i < Row_Number; ++i) {
			for (int j = 0; j < Col_Number; ++j) {
				 G._matrix[i][j] = G_._matrix[j][i];
				//G._matrix[i][j] = G_._matrix[j][i];
			}
		}
		//MatrixForm_to_Generator(G);// Polar Code不能做systematic form!
		break;
	case Polar_256_192_AStar_192_128:   //可隨意調整長度!
		Col_Number = Polar_Astar_Outer_N;
		Row_Number = Polar_Astar_Outer_K;
		Title = "(256,n) Polar Code concatenate (n,128) Random Code";
		CreateRandomCode(Row_Number, Col_Number);
		ReadFile_GeneratorMatrix("RandomCode.txt", G);
		MatrixForm_to_Generator(G);
		Inner_Col_Number = Polar_Astar_Inner_N;
		Inner_Row_Number = Col_Number;
		G_.Building_Empty_Matrix(Inner_Col_Number, Inner_Row_Number);
		Polar_Code_G_Generator(G_, frozen, non_frozen, Channel_Parameter);
		G_Inner.Building_Empty_Matrix(Inner_Row_Number, Inner_Col_Number);
		for (int i = 0; i < Inner_Row_Number; ++i) {
			for (int j = 0; j < Inner_Col_Number; ++j) {
				G_Inner._matrix[i][j] = G_._matrix[j][i];
			}
		}
		break;
	case LDPC_256_192_AStar_192_128:
		Col_Number = 192;
		Row_Number = 128;
		CreateRandomCode(Row_Number, Col_Number);
		ReadFile_GeneratorMatrix("RandomCode.txt", G);
		MatrixForm_to_Generator(G);
		Inner_Col_Number = 320;
		Inner_Row_Number = Col_Number;
		G_Inner.Col_number = 320;
		G_Inner.Row_number = 192;
		Title = "(256,192) LDPC Code concatenate (192,128) Random Code";
		ReadFile_H_Matrix_ver2("H_256_192.txt", H);
		Permutation_Seq.resize(Inner_Col_Number, 0);
		LDPC.H_G_convertor_G_P_I(H, Permutation_Seq, G_Inner);
		Punctured_Number = 64;
		for (int i = 0; i < G_Inner.Row_number; ++i) {    // Delete Col by Puncturing
			G_Inner._matrix.at(i).erase(G_Inner._matrix.at(i).begin(), G_Inner._matrix.at(i).begin() + (Punctured_Number));
		}
		G_Inner.Col_number -= Punctured_Number;
		Inner_Col_Number = 256;
		/*
		for (int i = 0; i < G.Row_number; ++i) {
			for (int j = 0; j < 64; ++j) {
				cout << (int)G._matrix[i][j];
			}
			cout << endl;
		}
		cout << endl;
		for (int i = 0; i < G.Row_number; ++i) {
			for (int j = 0; j < 192; ++j) {
				cout << (int)G._matrix[i][j + 64];
			}
			cout << endl;
		}
		cout << 
		system("pause");*/
		break;

	case Hamming_255_187_Astar_187_128:
		Col_Number = 187;              // 187
		Row_Number = 128;
		Short_Inner_Col_Number = 15;   // 15 
		Short_Inner_Row_Number = 11;   // 11 
		Title = "Outer RC(187, 128) Inner 17 Hamming Code(255, 187)";
		CreateRandomCode(Row_Number, Col_Number);
		ReadFile_GeneratorMatrix("RandomCode.txt", G);
		MatrixForm_to_Generator(G);
		HammingCode_Generator_Matrix(G_Inner,H, poly_4, Short_Inner_Col_Number);
		break;

	case RM_256_192_Astar_192_128:
		Col_Number = 176;              // 176
		Row_Number = 128;
		Short_Inner_Col_Number = 16;   // 16 
		Short_Inner_Row_Number = 11;   // 11 
		Title = "Outer RC (176,128) Inner 16 RM Code (256,176)";
		CreateRandomCode(Row_Number, Col_Number);
		ReadFile_GeneratorMatrix("RandomCode.txt", G);
		MatrixForm_to_Generator(G);
		RM_16_11_Generator_Matrix(G_Inner);
		break;

	case RM_256_168_Astar_168_128:
		Col_Number = 168;              // 176
		Row_Number = 128;
		Short_Inner_Col_Number = 64;   // 16 
		Short_Inner_Row_Number = 42;   // 11 
		Title = "Outer RC (168,128) Inner 4 RM Code (256,168)";
		CreateRandomCode(Row_Number, Col_Number);
		ReadFile_GeneratorMatrix("RandomCode.txt", G);
		MatrixForm_to_Generator(G);
		RM_order3_Generator_Matrix(6, G_Inner, G_);
		break;
	case Turbo_256_128_with_2LDPC:
		Title = "(256,128) Turbo Code with 2 (192,128) LDPC Code";
		// Candidates:
		// QC_LDPC_g10_H_192_64_hm_6_2.txt
		// QC_LDPC_g10_H_192_64_hm_6_2_ver2.txt
		// QC_LDPC_g10_H_192_64_hm_6_2_ver3.txt
		// QC_LDPC_192_128.txt
		// QC_LDPC_192_128(1,3).txt
		// H_64x192_z8.txt
		// QC_LDPC_H_192_64_hm_6_2.txt
		// QC_LDPC_H_192_64_hm_6_2_o.txt
		// QC_LDPC_H_192_64_hm_6_2_x.txt
		// QC_LDPC_H_192_64_hm_6_2_test.txt
		// QC_LDPC_H_189_63_hm_9_3.txt
		// QC_LDPC_H_192_64_hm_12_4.txt

		// girth 6 : 1: QC , 2: x
		// girth 10: 1: ver2, 2: ver3

		// Encoder 2
		ReadFile_GeneratorMatrix("QC_LDPC_g10_H_192_64_hm_6_2_ver2.txt", H);     
		// Encoder 1
		ReadFile_GeneratorMatrix("QC_LDPC_g10_H_192_64_hm_6_2_ver3.txt", H_);

		//Cycle4_Check(H);

		//Col_Number = H.Col_number;
		//Row_Number = H.Col_number - H.Row_number;
		Col_Number = 256;
		Row_Number = 128;

		Permutation_Seq.resize(Col_Number, 0);
		Permutation_Seq_.resize(Col_Number, 0);
		LDPC.fullrank = 1;
		LDPC.H_G_convertor_G_I_P_ver2(H, Permutation_Seq, G);
		LDPC.fullrank = 1;
		LDPC.H_G_convertor_G_I_P_ver2(H_, Permutation_Seq_, G_);
		break;
	case Turbo_256_128_with_2LDPC_Punc:
		Title = "(256,128) Turbo Code with 2 (192,128) LDPC Code, Punctured version";
		ReadFile_H_Matrix_ver2("H_256_128.txt", H);
		ReadFile_H_Matrix_ver2("H_256_128.txt", H_);
		Col_Number = 256;
		Row_Number = 128;
		G.Col_number = 256;
		G.Row_number = 128;
		LDPC.fullrank = 1;
		Permutation_Seq.resize(Col_Number, 0);
		LDPC.H_G_convertor_G_I_P_ver2(H, Permutation_Seq, G);
		Permutation_Seq_.resize(Col_Number, 0);
		LDPC.H_G_convertor_G_I_P_ver2(H_, Permutation_Seq, G_);
		//system("pause");
		//LDPC.H_G_convertor_G_P_I(H, Permutation_Seq, G);
		break;

	case Turbo_256_128_with_2LDPC_ver2:
		Title = "(256,128) Turbo Code with (160,128) LDPC Code & (224,128) LDPC Code";

		// QC_LDPC_g6_H_160_32_hm_10_2
		// QC_LDPC_g6_H_224_96_hm_7_3
		// QC_LDPC_H_224_96_hm_7_3
		// QC_LDPC_H_160_32_hm_10_2
		// Encoder 2
		ReadFile_GeneratorMatrix("QC_LDPC_g6_H_224_96_hm_7_3.txt", H);

		// Encoder 1
		ReadFile_GeneratorMatrix("QC_LDPC_g6_H_160_32_hm_10_2.txt", H_);

		//Col_Number = 256;
		//Row_Number = 128;
		Col_Number = 256;
		Row_Number = 128;

		Permutation_Seq.resize(H.Col_number, 0);  // Encoder 2
		Permutation_Seq_.resize(H_.Col_number, 0); // Encoder 1
		//cout << H.Col_number << "," << H.Row_number << "/" << H_.Col_number << "," << H_.Row_number << endl;
		LDPC.fullrank = TRUE;
		LDPC.H_G_convertor_G_I_P_ver2(H, Permutation_Seq, G);
		LDPC.fullrank = TRUE;
		LDPC.H_G_convertor_G_I_P_ver2(H_, Permutation_Seq_, G_);
		//system("pause");
		break;
	case Turbo_320_128_with_2LDPC:
		Title = "(320,128) Turbo Code with (224,128) LDPC Code & (224,128) LDPC Code";

		// Encoder 2
		ReadFile_GeneratorMatrix("QC_LDPC_H_224_96_hm_7_3.txt", H);

		// Encoder 1
		ReadFile_GeneratorMatrix("QC_LDPC_H_224_96_hm_7_3.txt", H_);

		Col_Number = 320;
		Row_Number = 128;
		//Col_Number = H.Col_number;
		//Row_Number = H.Col_number - H.Row_number;

		Permutation_Seq.resize(H.Col_number, 0);  // Encoder 2
		Permutation_Seq_.resize(H_.Col_number, 0); // Encoder 1
		//cout << H.Col_number << "," << H.Row_number << "/" << H_.Col_number << "," << H_.Row_number << endl;
		LDPC.fullrank = TRUE;
		LDPC.H_G_convertor_G_I_P_ver2(H, Permutation_Seq, G);
		LDPC.fullrank = TRUE;
		LDPC.H_G_convertor_G_I_P_ver2(H_, Permutation_Seq_, G_);
		//system("pause");
		break;
	case RandomCode_with_LDPC:
		Row_Number = 96;
		Col_Number = 192;
		CreateRandomCode(Row_Number, Col_Number);
		ReadFile_GeneratorMatrix("RandomCode.txt", G);
		MatrixForm_to_Generator(G);
		Convert_SystematicG_to_H(H, G);
		break;
	case Turbo_256_128_RC_LDPC_192_128:
		Title = "(256,128) Turbo Code with RC & LDPC Code (192,128)";
		// Encoder 1
		ReadFile_GeneratorMatrix("QC_LDPC_192_128.txt", H_);

		Col_Number = 256;
		Row_Number = 128;
		Permutation_Seq_.resize(H_.Col_number, 0);
		LDPC.H_G_convertor_G_I_P_ver2(H_, Permutation_Seq_, G_);
		// Encoder 2
		CreateRandomCode(Row_Number, 192);
		ReadFile_GeneratorMatrix("RandomCode.txt", G);
		MatrixForm_to_Generator(G);
		Convert_SystematicG_to_H(H, G);
		break;

	case Polar_with_Astar_middle_index:
		Title = "Turbo Code with Astar in the middle";
		Col_Number = 1024;
		Row_Number = 512;
		G_.Building_Empty_Matrix(Col_Number, Row_Number + Outer_Astar_Length / 2);
		Polar_Code_G_Generator(G_, frozen, non_frozen, Channel_Parameter);
		G.Building_Empty_Matrix(Row_Number + Outer_Astar_Length / 2, Col_Number);
		for (int i = 0; i < Row_Number + Outer_Astar_Length / 2; ++i) {
			for (int j = 0; j < Col_Number; ++j) {
				G._matrix[i][j] = G_._matrix[j][i];
				//G._matrix[i][j] = G_._matrix[j][i];
			}
		}
		CreateRandomCode(Outer_Astar_Length / 2, Outer_Astar_Length);
		ReadFile_GeneratorMatrix("RandomCode.txt", G_Inner);
		MatrixForm_to_Generator(G_Inner);
		//Col_Number = 256;
		break;

	default:
		cout << " Please check your choice !" << endl;
		system("pause");
		exit(1);
		break;
	}
	if (Col_Number < 257) {
		cout
			<< "\n " << Title << "\n"
			<< " Generator Matrix in Systematic Form: \n"
			<< endl;
		G.Show_matrix(1);
		system("pause");
	}
	if (Inner_Col_Number != 0) Code_Rate = (double)Row_Number / ((double)Inner_Col_Number);
	//if(CRC_Check==TRUE) Code_Rate = ((double)(Row_Number-CRC_length)) / ((double)(Col_Number));
	else if (Short_Inner_Col_Number != 0)  Code_Rate = (double)Row_Number*(double)Short_Inner_Row_Number / ((double)Col_Number) / ((double)Short_Inner_Col_Number);
	else Code_Rate = ((double)Row_Number) / ((double)Col_Number);
	//Code_Rate = ((double)Row_Number) / ((double)Col_Number);
	//cout << Code_Rate << "," << Row_Number << "," << Inner_Col_Number << endl;
	//system("pause");
}

void CODE::CRC_Desicion(vector<__int8>& CRC) {
	if (CRC_length == 16) CRC = CRC_16;
	else if (CRC_length == 8) CRC = CRC_8;  //  { 1,0,0,1,1,0,0,0,1 } , bluetooth: { 1,1,0,1,0,0,1,1,1 }
	else if (CRC_length == 10) CRC = CRC_10;
	else if (CRC_length == 12) CRC = CRC_12;
	else if (CRC_length == 14) CRC = CRC_14;
	else if (CRC_length == 24) CRC = CRC_24;
	else cout << "Wrong CRC length!" << endl;
}

vector<__int8> CODE::CRC_Generator(vector<__int8> Message, vector<__int8> variant) {
	int digit = 0;
	int Boundary = Message.size() - 1;
	vector<__int8> Dividend = Message;
    int Variant_Size = variant.size();
	//cout << Variant_Size << "/" << Boundary << endl;

	Dividend.insert(Dividend.end(), Variant_Size - 1, 0);

	while (digit <= Boundary) {
		if (Dividend.at(digit) == 1) {
			for (int i = 0; i < Variant_Size; ++i) {
				Dividend.at(i + digit) ^= variant.at(i);
			}
		}
		//for (int j = 0; j < Dividend.size(); ++j) std::cout << (int)Dividend.at(j) << " ";
		//std::cout << std::endl;
		++digit;
	}
	//cout << digit << endl;
	vector<__int8> Return_CRC;
	Return_CRC.assign(Dividend.begin() + Boundary + 1, Dividend.end());
	return Return_CRC;
}


bool CODE::CRC_Examination(vector<__int8> Message_and_CRC, vector<__int8> variant) {
	__int8 digit = 0;
	__int8 Boundary = Message_and_CRC.size() - variant.size();
	vector<__int8> Dividend = Message_and_CRC;
	__int8 Variant_Size = variant.size();

	while (digit <= Boundary) {
		if (Dividend.at(digit) == 1) {
			for (int i = 0; i < Variant_Size; ++i) {
				Dividend.at(i + digit) ^= variant.at(i);
			}
		}
		//for (int j = 0; j < Dividend.size(); ++j) std::cout << (int)Dividend.at(j) << " ";
		//std::cout << std::endl;
		++digit;
	}
	//std::cout << "AAAA";
	--digit;
	bool Result = 1;
	while (digit < Message_and_CRC.size()) {
		if (Dividend.at(digit) == 1) {
			Result = 0;
			break;
		}
		//for (int j = 0; j < Dividend.size(); ++j) std::cout << (int)Dividend.at(j) << " ";
		//std::cout << std::endl;
		++digit;
	}
	return Result;
}


void MatrixForm_to_Generator(MATRIX<__int8> &G)
{
	GJ_Elimination(G);
	// 將做完高斯消去的 generator matrix 調整為 [ I | P ] 的矩陣型式
	size_t Register(0), col_number(0);
	for (size_t i(0); i < G.Row_number; ++i){	//檢查前半部是否為 identity matrix
		if (G._matrix[i][i] != 1){
			col_number = 0;
			//同一列，由下一行開始，碰到第一個此列為 "1" 的行就交換 column index
			for (size_t j(i + 1); j < G.Col_number; ++j){
				if (G._matrix[i][j] == 1){
					col_number = j;
					break;
				}
			}
			//Exchange those columns
			Register = 0;
			for (size_t k(0); k < G.Row_number; ++k){
				Register = G._matrix[k][i];
				G._matrix[k][i] = G._matrix[k][col_number];
				G._matrix[k][col_number] = Register;
			}
		}
	}
};

//By 張家輔學長
void BCH_Generator_Matrix(size_t G1[], MATRIX<__int8> &G)
{
	size_t
		reg(0),
		Col_Number(G.Col_number),
		Row_Number(G.Row_number);

	/* 填入 G 第一列G1，第一個數字作邏輯右移，找到第一個 "1" 的位置 ，用 k 儲存 */
	for (size_t i(0); i < 30; ++i){
		if (((G1[i / 30] >> (29 - (i % 30))) & 1) == 1){
			reg = i;
			break;
		}
	}
	for (size_t i(reg); i <= Col_Number - Row_Number + reg - 1; ++i)
		G._matrix[0][i - reg] = (G1[i / 30] >> (29 - (i % 30))) & 1;

	/* 將第一列右移，當作下一列，重複至將 G 的所有列都填完，再作 Gaussian elimination  */
	/* 得到 G
	G =
	｢   |   | 0 |
	|   |   | . |
	| I | P | . |
	|   |   | . |
	|   |   | 0 ｣
	*/
	for (size_t i(0); i < Col_Number - Row_Number; ++i){
		if (G._matrix[0][i] == 1)
			for (size_t j(1); j < Row_Number; ++j)
				G._matrix[j][j + i] = 1;
	}
	GJ_Elimination(G);
}

void Extended_BCH_Generator_Matrix(size_t G1[], MATRIX<__int8> &G)
{
	size_t
		reg(0),
		Col_Number(G.Col_number),
		Row_Number(G.Row_number);

	/* 填入 G 第一列G1，第一個數字作邏輯右移，找到第一個 "1" 的位置 ，用 k 儲存 */
	for (size_t i(0); i < 30; ++i){
		if (((G1[i / 30] >> (29 - (i % 30))) & 1) == 1){
			reg = i;
			//cout << "reg = " << reg << endl;
			break;
		}
	}
	for (size_t i(reg); i <= Col_Number - Row_Number + reg - 1; ++i)
		G._matrix[0][i - reg] = (G1[i / 30] >> (29 - (i % 30))) & 1;

	/* 將第一列右移，當作下一列，重複至將 G 的所有列都填完，再作 Gaussian elimination  */
	/* 得到 G
	G =
	[   |   | 0 ]
	[   |   | . ]
	[ I | P | . ]
	[   |   | . ]
	[   |   | 0 ]
	*/
	for (size_t i(0); i < Col_Number - Row_Number; ++i){
		if (G._matrix[0][i] == 1)
			for (size_t j(1); j < Row_Number; ++j)
				G._matrix[j][j + i] = 1;
	}
	GJ_Elimination(G);
	
	/* 新增 Temp_H */
	/*
	Temp_H 為 (N - K1) * N matrix，其第一列與第一行均為 0 向量，
	將 P 的 transpose 放置在以 (2,2) 為起始位置的部分，其餘填零
	Temp_H =
	｢ 0 |  0 .... 0  |
	| . |            |
	| . |            |
	| . |    P^T     |
	| . |            |
	| 0 |            ｣
	*/
	MATRIX<__int8> Temp_H(Col_Number - Row_Number, Col_Number);
	for (size_t i(1); i <= Row_Number; ++i)
		for (size_t j(1); j < Col_Number - Row_Number; ++j)
			Temp_H._matrix[j][i] = G._matrix[i - 1][j + Row_Number - 1];
	/* 將 Temp_H 第一列全填 1，後半部的零矩陣，填入identity matrix，再作 Gaussian elimination */
	for (size_t i(0); i < Col_Number; ++i)	Temp_H._matrix[0][i] = 1; // 第一列全填 1
	for (size_t i(Row_Number + 1); i < Col_Number; ++i)	Temp_H._matrix[i - Row_Number][i] = 1; // 填入identity matrix
	GJ_Elimination(Temp_H);
	/* 將 Temp_H 後半部轉置再填入 Generator matrix 的後半部 */
	for (size_t i(0); i < Row_Number; ++i)
		for (size_t j(0); j < Col_Number - Row_Number; ++j)
			G._matrix[i][j + Row_Number] = Temp_H._matrix[j][i + Col_Number - Row_Number];
}

void Systematic_Linear_Block_Code_Encoder(MATRIX<__int8> &G, vector<__int8> &message_seq, vector<__int8> &output_codeword_seq)
{// systematic only
	vector<__int8> codeword_seq(G.Col_number, 0);

	// computation for parity part
	/*for (size_t i(0); i < G.Row_number; ++i){
		if (message_seq.at(i) == 1) {
			codeword_seq.at(i) = 1;
			for (size_t j(G.Row_number); j < G.Col_number; ++j) {
				if (G._matrix[i][j] == 1)
					codeword_seq.at(j) ^= 1;
			}
		}
	}*/

	//2019
	/*for (size_t index(0); index < G.Row_number; ++index) {
		if (message_seq.at(index) == 1) {
			auto output = (G._matrix[index] ^ codeword_seq);
			for (size_t i(0); i < G.Col_number; ++i) {
				codeword_seq.at(i) = output[i];
			}
		}
	}*/

	//2019/07/10 new method，this method is more efficient and simpler.
	for (size_t index(0); index < G.Row_number; ++index) {
		if (message_seq.at(index) == 1) {
			transform(
				codeword_seq.begin(),
				codeword_seq.end(),
				G._matrix[index].begin(),
				codeword_seq.begin(),
				std::bit_xor<__int8>());
		}
	}

	output_codeword_seq = codeword_seq;
} // end Linear_Block_Code_Encoder()

void Show_All_Codeword(MATRIX<__int8> G)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		message_state_number((unsigned int)pow(2, message_length));
	vector <__int8>
		message_seq(message_length, 0),
		codeword_seq(codeword_length, 0);

	//Show All Code Word
	std::cout << "\n Code Word：";
	for (size_t i(0); i < message_state_number; ++i){
		//給定 message sequence
		for (size_t j(0); j < message_length; ++j)
			message_seq.at(message_length - j - 1) = (i >> j) & 1;

		//初始化 codeword sequence
		codeword_seq.clear();
		Systematic_Linear_Block_Code_Encoder(G, message_seq, codeword_seq);
		std::cout << std::endl;
		for (size_t i(0); i < message_length; ++i) std::cout << " " << message_seq.at(i);

		std::cout << " >> ";
		for (size_t i(0); i < codeword_length; ++i) std::cout << " " << codeword_seq.at(i);
	};
	//Show Generator matrix
	cout << "\n Generator Matrix：" << endl;
	for (size_t i(0); i < G.Row_number; ++i){
		for (size_t j(0); j < G.Col_number; ++j)
			cout << " " << G._matrix[i][j];
		cout << endl;
	}
	cout << endl;
}//end Show_All_Codeword() 沒用到

void ReadFile_GeneratorMatrix(string name,MATRIX<__int8> &G)
{
	fstream fp;
	fp.open(name, ios::out | ios::in);
	if (!fp) cout << "Fail to open file: " << endl;
	char ch;
	size_t i(0), j(0), Row_number(0), Col_number(0);
	fp >> Row_number;
	fp >> Col_number;
	G.Building_Empty_Matrix(Row_number, Col_number);

	fp.get(ch);
	while (fp.get(ch)) {
		if (ch == '\n') {
			++i;
			j = 0;
		}
		else if (ch == ' '); // Chou yin
		else {
			G._matrix[i][j] = ch - '0';
			++j;
		}
	}
	fp.close();
};

void CreateRandomCode(size_t message_length, size_t codeword_length) {
	fstream fp;
	fp.open("RandomCode.txt", ios::out | ios::trunc);
	if (!fp) cout << "Fail to open file: " << endl;
	fp << message_length << "\n";
	fp << codeword_length << "\n";

	std::default_random_engine generator;
	std::bernoulli_distribution distribution(0.5);

	for (size_t i(0); i < message_length; ++i) {
		for (size_t j(0); j < codeword_length; ++j) {
			fp << distribution(generator);
		}
		if (i == (message_length - 1));
		else fp << "\n";
	}
	fp.close();
}

void CreateRandomCode(size_t row_start_position, size_t row_end_position, double probability, MATRIX<__int8> &H) {
	std::default_random_engine generator;
	std::bernoulli_distribution distribution(probability);
	for (size_t i(row_start_position); i < row_end_position; ++i) {
		for (size_t j(0); j < H.Col_number; ++j) {
			H._matrix[i][j] = distribution(generator);
		}
	}
}

void FillRandomBits(size_t row_start_position, MATRIX<__int8> &H) {
	std::default_random_engine generator;
	std::bernoulli_distribution distribution(0.5);
	for (size_t i(row_start_position); i < H.Row_number; i++) {
		for (size_t j(0); j < H.Col_number; j++) {
			H._matrix[i][j] = distribution(generator);
		}
	}
}

void GenerateParityMatrix(MATRIX<__int8> &G) {
	MATRIX<__int8> P_matrix;
	P_matrix.Building_Empty_Matrix(G.Col_number - G.Row_number, G.Col_number);

	//identity part of Parity matrix
	for (size_t i(0); i < G.Col_number - G.Row_number; ++i) {
		P_matrix._matrix[i][i + G.Row_number] = 1;
	}
	for (size_t i(0); i < G.Row_number; ++i) {
		for (size_t j(0); j < G.Col_number-G.Row_number; ++j) {
			P_matrix._matrix[i][j] = G._matrix[j][G.Col_number - G.Row_number + i];
		}
	}

	fstream fp;
	fp.open("Polar_ParityMatrix.txt", ios::out | ios::trunc);
	if (!fp) cout << "Fail to open file: " << endl;
	fp << G.Row_number << "\n";
	fp << G.Col_number << "\n";

	for (size_t i(0); i < G.Col_number-G.Row_number; ++i) {
		for (size_t j(0); j < G.Col_number; ++j) {
			fp << P_matrix._matrix[i][j];
			if (j == (G.Col_number - 1)) fp << "\n";
		}
	}
	fp.close();
}

void HammingCode_Generator_Matrix(MATRIX<__int8> &G, MATRIX<__int8> &H,vector<__int8> &generator_polynomial, __int8 codelength) {
	__int8 parity_check_bit_length = generator_polynomial.size() - 1;
	__int8 row_number = codelength - parity_check_bit_length;
	G.Building_Empty_Matrix(row_number, codelength);
	vector<__int8> Dividend(codelength, 0);

	for (int row = 0; row < row_number; ++row) {
		G._matrix[row][row] = 1;
		Dividend.at(row) = 1;
		for (int shift_register = 0; shift_register < row_number; shift_register++) {
			if (Dividend.at(shift_register) == 1) {
				for (int j = 0; j < generator_polynomial.size(); ++j)  Dividend.at(shift_register + j) ^= generator_polynomial.at(j);
			}
		}
		for (int j = row_number; j < codelength; ++j) {
			G._matrix[row][j] = Dividend.at(j);
			Dividend.at(j) = 0;
		}
	}

	H.Building_Empty_Matrix(codelength - row_number, codelength);
	for (int row = 0; row < codelength - row_number; ++row) {
		H._matrix[row][row] = 1;
		for (int col = 0; col < row_number; ++col) {
			H._matrix[row][col + parity_check_bit_length] = G._matrix[col][row + row_number];
		}
	}
	/*
	for (int i = 0; i < G.Row_number; ++i) {
		for (int j = 0; j < G.Col_number; ++j) {
			cout << (int)G._matrix[i][j];
		}
		cout << endl;
	}
	cout << endl;
	cout << endl;
	for (int i = 0; i < H.Row_number; ++i) {
		for (int j = 0; j < H.Col_number; ++j) {
			cout << (int)H._matrix[i][j];
		}
		cout << endl;
	}
	system("pause");*/

}

void RM_16_11_Generator_Matrix(MATRIX<__int8> &G) {
	G.Building_Empty_Matrix(11, 16);
	G._matrix.at(0) = RM_16_v0;
	G._matrix.at(1) = RM_16_v4;
	G._matrix.at(2) = RM_16_v3;
	G._matrix.at(3) = RM_16_v2;
	G._matrix.at(4) = RM_16_v1;
	G._matrix.at(5) = RM_16_v3v4;
	G._matrix.at(6) = RM_16_v2v4;
	G._matrix.at(7) = RM_16_v1v4;
	G._matrix.at(8) = RM_16_v2v3;
	G._matrix.at(9) = RM_16_v1v3;
	G._matrix.at(10) = RM_16_v1v2;
	/*
	for (int i = 0; i < G.Row_number; ++i) {
		for (int j = 0; j < G.Col_number; ++j) {
			cout << (int)G._matrix[i][j];
		}
		cout << endl;
	}
	cout << endl;
	system("pause");*/
}

void RM_32_6_Generator_Matrix(MATRIX<__int8> &G) {
	G.Building_Empty_Matrix(6, 32);
	G._matrix.at(0) = RM_32_v0;
	G._matrix.at(1) = RM_32_v5;
	G._matrix.at(2) = RM_32_v4;
	G._matrix.at(3) = RM_32_v3;
	G._matrix.at(4) = RM_32_v2;
	G._matrix.at(5) = RM_32_v1;
	
}

void Reed_Muller_Code_Generator_Matrix(MATRIX<__int8> &G, __int8 r, __int8 m) {   // Generate RM(r,m) -> (n,k) RM Code , n = 2^m, k = 1 + (m 1) + ... + (m r) 
	int n = pow(2, m);

	G.Building_Empty_Matrix(6, 32);
	

}


void RM_order3_Generator_Matrix(int m, MATRIX<__int8> &G, MATRIX<__int8> &G_) {
	// r = 3
	int r = 3;         // order is 3
	int n = pow(2, m); // length of codeword
	int k = 0;         // initialize length of message 
	vector<__int8> order_length_record;

	for (int i = r; i > -1; --i) {
		int temp = combination(m, r - i);
		k += temp;
		order_length_record.push_back(temp);
	}
	G.Building_Empty_Matrix(k, n);
	//G_.Building_Empty_Matrix(k, n);
	G_.Building_Empty_Matrix(k, 0);
	// r = 0
	vector<__int8> All_one_vector(n, 1);
	G._matrix.at(0) = All_one_vector;

	// 若長度調整, 這個地方要改!!!
	MATRIX<__int8> Order1_vector, Order2_vector, Order3_vector, Order1_set, Order2_set, Order3_set;
	Order1_vector.Building_Empty_Matrix(m, n);
	Order1_vector._matrix.at(0) = RM_64_v1;
	Order1_vector._matrix.at(1) = RM_64_v2;
	Order1_vector._matrix.at(2) = RM_64_v3;
	Order1_vector._matrix.at(3) = RM_64_v4;
	Order1_vector._matrix.at(4) = RM_64_v5;
	Order1_vector._matrix.at(5) = RM_64_v6;

	Order2_vector.Building_Empty_Matrix(order_length_record.at(2), n);
	Order3_vector.Building_Empty_Matrix(order_length_record.at(3), n);
	Order1_set.Building_Empty_Matrix(order_length_record.at(1), n);
	Order2_set.Building_Empty_Matrix(order_length_record.at(2), n);
	Order3_set.Building_Empty_Matrix(order_length_record.at(3), n);
	// end
	vector<string> R0, R1, R2, R3;
	string StringTemp;
	// r = 1
	int NumberOfOrder2 = 0, NumberOfOrder3 = 0;
	for (int order1 = 0; order1 < m; ++order1) {
		Order1_set._matrix.at(order1) = DefineSets_order1(order1, m);
		
		//string StringTest = "(v" + to_string(1) + ", v" + to_string(2) + ")";
		StringTemp = "(v" + to_string(order1+1) + ")";
		R1.push_back(StringTemp);
		// r = 2
		for (int order2 = order1 + 1; order2 < m; ++order2) {
			Order2_vector._matrix.at(NumberOfOrder2) = Vector_Mulitiplication(Order1_vector._matrix.at(order1), Order1_vector._matrix.at(order2));
			Order2_set._matrix.at(NumberOfOrder2++) = DefineSets_order2(order1, order2, m);

			StringTemp = "(v" + to_string(order1 + 1) + ",v" + to_string(order2 + 1) + ")";
			R2.push_back(StringTemp);
			// r = 3
			for (int order3 = order2 + 1; order3 < m; ++order3) {
				Order3_vector._matrix.at(NumberOfOrder3) = 
					Vector_Mulitiplication(Order1_vector._matrix.at(order1), Order1_vector._matrix.at(order2), Order1_vector._matrix.at(order3));
				Order3_set._matrix.at(NumberOfOrder3++) = DefineSets_order3(order1, order2, order3, m);

				StringTemp = "(v" + to_string(order1 + 1) + ",v" + to_string(order2 + 1) + ",v" + to_string(order3 + 1) + ")";
				R3.push_back(StringTemp);
			}
		}
	}

	// Construct Generator matrix G
	// r = 0
	int GCounter = 1, G_Counter = 1;
	// r = 1
	for (int count = 0; count < order_length_record.at(1); count++) {
		G._matrix.at(GCounter++) = Order1_vector._matrix.at(count);
		G_._matrix.at(G_Counter++) = Order1_set._matrix.at(count);
	}
	// r = 2
	for (int count = 0; count < order_length_record.at(2); count++) {
		G._matrix.at(GCounter++) = Order2_vector._matrix.at(count);
		G_._matrix.at(G_Counter++) = Order2_set._matrix.at(count);
	}
	// r = 3
	for (int count = 0; count < order_length_record.at(3); count++) {
		G._matrix.at(GCounter++) = Order3_vector._matrix.at(count);
		G_._matrix.at(G_Counter++) = Order3_set._matrix.at(count);
	}
	/*
	for (int i = 0; i < 42; ++i) {
		cout << G_._matrix.at(i).size() << " ";
	}
	cout << endl;
	*/
	// Add the affecting rows into Independent Determinations
	// r = 0 
	__int8 Row_Check_Start = order_length_record.at(0);
	for (int set = 0; set < n; set++) {
		G_._matrix.at(0).push_back(set);
		for (int row_check = Row_Check_Start; row_check < k; row_check++) {
			if (G._matrix[row_check][set] == 1) G_._matrix.at(0).push_back(row_check*(-1));
		}
		G_._matrix.at(0).push_back(CHAR_MAX);
	}

	// r = 1 
	Row_Check_Start = order_length_record.at(0) + order_length_record.at(1);
	vector<__int8> Set_Temp;
	int check_temp;

	for (int row = order_length_record.at(0); row < k; row++) {
		if (row == order_length_record.at(0) + order_length_record.at(1))
			Row_Check_Start = order_length_record.at(0) + order_length_record.at(1) + order_length_record.at(2);
		if (row == order_length_record.at(0) + order_length_record.at(1) + order_length_record.at(2))
			Row_Check_Start = k;
		int traversal = 0;
		while (traversal < G_._matrix.at(row).size()) {
			while (G_._matrix[row][traversal] != CHAR_MAX) {
				Set_Temp.push_back(G_._matrix[row][traversal]);
				++traversal;
			}
			for (int check_row = Row_Check_Start; check_row < k; check_row++) {
				check_temp = 0;
				for (int set_traversal = 0; set_traversal < Set_Temp.size(); set_traversal++) {
					check_temp ^= G._matrix[check_row][Set_Temp.at(set_traversal)];
				}
				if (check_temp == 1) {
					G_._matrix.at(row).insert(G_._matrix.at(row).begin() + traversal, check_row*(-1));
					++traversal;
				}
			}
			++traversal;
			Set_Temp.clear();
		}
	}

	// Delete Overlapped Calculations
	
	for (int row = 0; row < k; ++row) {
		int Set_Start = 0, accumulated_negative_number = 0;
		for (int Set_traversal = 0; Set_traversal < G_._matrix.at(row).size(); ++Set_traversal) {
			if (G_._matrix[row][Set_traversal] < 0) accumulated_negative_number++;
			if (accumulated_negative_number > 1) {
				int Interval = 0;
				while (G_._matrix[row][Set_Start + Interval] != CHAR_MAX) {
					Interval++;
				}
				G_._matrix.at(row).erase(G_._matrix.at(row).begin() + Set_Start, G_._matrix.at(row).begin() + Set_Start + Interval+ 1);
				Set_traversal = Set_Start;
				accumulated_negative_number = 0;
			}
			if (G_._matrix[row][Set_traversal] == CHAR_MAX) {
				Set_Start = Set_traversal + 1;
				accumulated_negative_number = 0;
			}
		}
	}


	/*
	for (int i = 0; i < 42; ++i) {
		cout << "row " << i << " : " << endl;
		for (int j = 0; j < G_._matrix.at(i).size(); ++j) {
			if (G_._matrix[i][j] == CHAR_MAX) cout << " / ";
			else if (G_._matrix[i][j] < 0) cout << "row" << G_._matrix[i][j] * (-1) <<" ";
			else cout << (int)G_._matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	system("pause");
	*/
	// Start Printing whatever
	/*
	StringTemp = "(v0)";
	R0.push_back(StringTemp);
	R0.insert(R0.end(), R1.begin(), R1.end());
	R0.insert(R0.end(), R2.begin(), R2.end());
	R0.insert(R0.end(), R3.begin(), R3.end());

	string filename("RM_64_42.txt");
	ClearFile(filename);
	fstream fp;
	fp.open(filename, ios::out | ios::app);
	for (int row = 0; row < k; row++) {
		if (row == 0) fp << "\nOrder 0: " << "\n\n";
		else if (row == 1) fp << "\nOrder 1: " << "\n\n";
		else if (row == 7) fp << "\nOrder 2: " << "\n\n";
		else if (row == 22) fp << "\nOrder 3: " << "\n\n";
		fp << "Row" << row << " " << R0.at(row) << ":\n";
		int SetCount = 1;
		fp << "Set1 ";
		for (int j = 0; j < G_._matrix.at(row).size(); ++j) {
			if (G_._matrix[row][j] == CHAR_MAX) {
				fp << "\n";
				if (j != G_._matrix.at(row).size() - 1) fp << "Set" << ++SetCount << " ";
			}
			else if (G_._matrix[row][j] < 0) fp << "LLR" << (int)G_._matrix[row][j] * (-1) << " ";
			else fp << (int)G_._matrix[row][j] << " ";
		}
		fp << "\n";
	}

	fp.close();
	system("pause");
	*/
}

double combination(double n, double k) {
	if (k == 0) return 1;
	else return (n / k)*combination(n - 1, k - 1);
}

vector<__int8> Vector_Mulitiplication(vector<__int8> vec1, vector<__int8> vec2) {
	vector<__int8> Return_Vec(vec1.size(), 0);
	for (int i = 0; i < vec1.size(); ++i) {
		if (vec1.at(i) == 1 && vec2.at(i) == 1) Return_Vec.at(i) = 1;
	}
	return Return_Vec;
}

vector<__int8> Vector_Mulitiplication(vector<__int8> vec1, vector<__int8> vec2, vector<__int8> vec3) {
	vector<__int8> Return_Vec(vec1.size(), 0);
	for (int i = 0; i < vec1.size(); ++i) {
		if ((vec1.at(i) == 1 && vec2.at(i) == 1) && vec3.at(i) == 1) Return_Vec.at(i) = 1;
	}
	return Return_Vec;
}

vector<__int8> Return_S_order1(int i1) {
	vector<__int8> S;
	i1 = pow(2, i1);
	for (int i = 0; i < 2; ++i) {
		S.push_back(i1*i);
	}
	return S;
}

vector<__int8> Return_E_order1(int i1, int m) {
	vector<__int8> E;
	for (int i = 0; i < m; ++i) {
		if (i == i1) continue;
		E.push_back(i);
	}
	return E;
}

vector<__int8> Return_Sc_order1(vector<__int8> E) {
	vector<__int8> Sc;
	// 若長度調整, 這個地方要改!!!
	for (int i = 0; i < E.size(); ++i) {
		E.at(i) = pow(2, E.at(i));
	}
	Sc.push_back(0);
	Sc.push_back(E.at(0));
	Sc.push_back(E.at(1));
	Sc.push_back(E.at(2));
	Sc.push_back(E.at(3));
	Sc.push_back(E.at(4));
	Sc.push_back(E.at(0) + E.at(1));
	Sc.push_back(E.at(0) + E.at(2));
	Sc.push_back(E.at(0) + E.at(3));
	Sc.push_back(E.at(0) + E.at(4));
	Sc.push_back(E.at(1) + E.at(2));
	Sc.push_back(E.at(1) + E.at(3));
	Sc.push_back(E.at(1) + E.at(4));
	Sc.push_back(E.at(2) + E.at(3));
	Sc.push_back(E.at(2) + E.at(4));
	Sc.push_back(E.at(3) + E.at(4));
	Sc.push_back(E.at(0) + E.at(1) + E.at(2));
	Sc.push_back(E.at(0) + E.at(1) + E.at(3));
	Sc.push_back(E.at(0) + E.at(1) + E.at(4));
	Sc.push_back(E.at(0) + E.at(2) + E.at(3));
	Sc.push_back(E.at(0) + E.at(2) + E.at(4));
	Sc.push_back(E.at(0) + E.at(3) + E.at(4));
	Sc.push_back(E.at(1) + E.at(2) + E.at(3));
	Sc.push_back(E.at(1) + E.at(2) + E.at(4));
	Sc.push_back(E.at(1) + E.at(3) + E.at(4));
	Sc.push_back(E.at(2) + E.at(3) + E.at(4));
	Sc.push_back(E.at(0) + E.at(1) + E.at(2) + E.at(3));
	Sc.push_back(E.at(0) + E.at(1) + E.at(2) + E.at(4));
	Sc.push_back(E.at(0) + E.at(1) + E.at(3) + E.at(4));
	Sc.push_back(E.at(0) + E.at(2) + E.at(3) + E.at(4));
	Sc.push_back(E.at(1) + E.at(2) + E.at(3) + E.at(4));
	Sc.push_back(E.at(0) + E.at(1) + E.at(2) + E.at(3) + E.at(4));
	sort(Sc.begin(), Sc.end());
	// end
	return Sc;
}

vector<__int8> DefineSets_order1(int i1, int m) {
	int n = pow(2, m);
	vector<__int8> S = Return_S_order1(i1);
	vector<__int8> Sc = Return_Sc_order1(Return_E_order1(i1, m));
	//vector<__int8> Return_Vector(n, 0);
	vector<__int8> Return_Vector; // , Set_Keeper;
	int NumberOfSets = n / S.size();
	for (int Set = 0; Set < NumberOfSets; Set++) {
		for (int index = 0; index < S.size(); ++index) {
			//Return_Vector.at(Set*S.size() + index) = S.at(index) + Sc.at(Set);
			Return_Vector.push_back(S.at(index) + Sc.at(Set));
		}
		Return_Vector.push_back(CHAR_MAX);
	}
	return Return_Vector;
}

vector<__int8> Return_S_order2(int i1, int i2) {
	vector<__int8> S;
	i1 = pow(2, i1);
	i2 = pow(2, i2);
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			S.push_back(i1*i + i2 * j);
		}
	}
	sort(S.begin(), S.end());
	return S;
}

vector<__int8> Return_E_order2(int i1, int i2, int m) {
	vector<__int8> E;
	for (int i = 0; i < m; ++i) {
		if (i == i1 || i == i2) continue;
		E.push_back(i);
	}
	return E;
}

vector<__int8> Return_Sc_order2(vector<__int8> E) {
	vector<__int8> Sc;
	// 若長度調整, 這個地方要改!!!
	for (int i = 0; i < E.size(); ++i) {
		E.at(i) = pow(2, E.at(i));
	}
	Sc.push_back(0);
	Sc.push_back(E.at(0));
	Sc.push_back(E.at(1));
	Sc.push_back(E.at(2));
	Sc.push_back(E.at(3));
	Sc.push_back(E.at(0) + E.at(1));
	Sc.push_back(E.at(0) + E.at(2));
	Sc.push_back(E.at(0) + E.at(3));
	Sc.push_back(E.at(1) + E.at(2));
	Sc.push_back(E.at(1) + E.at(3));
	Sc.push_back(E.at(2) + E.at(3));
	Sc.push_back(E.at(0) + E.at(1) + E.at(2));
	Sc.push_back(E.at(0) + E.at(1) + E.at(3));
	Sc.push_back(E.at(0) + E.at(2) + E.at(3));
	Sc.push_back(E.at(1) + E.at(2) + E.at(3));
	Sc.push_back(E.at(0) + E.at(1) + E.at(2) + E.at(3));
	sort(Sc.begin(), Sc.end());
	// end
	return Sc;
}

vector<__int8> DefineSets_order2(int i1, int i2, int m) {
	int n = pow(2, m);
	vector<__int8> S = Return_S_order2(i1, i2);
	vector<__int8> Sc = Return_Sc_order2(Return_E_order2(i1, i2, m));
	//vector<__int8> Return_Vector(n, 0);
	vector<__int8> Return_Vector;
	int NumberOfSets = n / S.size();
	for (int Set = 0; Set < NumberOfSets; Set++) {
		for (int index = 0; index < S.size(); ++index) {
			//Return_Vector.at(Set*S.size() + index) = S.at(index) + Sc.at(Set);
			Return_Vector.push_back(S.at(index) + Sc.at(Set));
		}
		Return_Vector.push_back(CHAR_MAX);
	}
	return Return_Vector;
}

vector<__int8> Return_S_order3(int i1, int i2, int i3) {
	vector<__int8> S;
	i1 = pow(2, i1);
	i2 = pow(2, i2);
	i3 = pow(2, i3);
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				S.push_back(i1 * i + i2 * j + i3 * k);
			}
		}
	}
	sort(S.begin(), S.end());
	return S;
}

vector<__int8> Return_E_order3(int i1, int i2, int i3, int m) {
	vector<__int8> E;
	for (int i = 0; i < m; ++i) {
		if (i == i1 || i == i2 || i == i3) continue;
		E.push_back(i);
	}
	return E;
}

vector<__int8> Return_Sc_order3(vector<__int8> E) {
	vector<__int8> Sc;
	// 若長度調整, 這個地方要改!!!
	for (int i = 0; i < E.size(); ++i) {
		E.at(i) = pow(2, E.at(i));
	}
	Sc.push_back(0);
	Sc.push_back(E.at(0));
	Sc.push_back(E.at(1));
	Sc.push_back(E.at(2));
	Sc.push_back(E.at(0) + E.at(1));
	Sc.push_back(E.at(0) + E.at(2));
	Sc.push_back(E.at(1) + E.at(2));
	Sc.push_back(E.at(0) + E.at(1) + E.at(2));
	sort(Sc.begin(), Sc.end());
	// end
	return Sc;
}

vector<__int8> DefineSets_order3(int i1, int i2, int i3, int m) {
	int n = pow(2, m);
	vector<__int8> S = Return_S_order3(i1, i2, i3);
	vector<__int8> Sc = Return_Sc_order3(Return_E_order3(i1, i2, i3, m));
	//vector<__int8> Return_Vector(n, 0);
	vector<__int8> Return_Vector;
	int NumberOfSets = n / S.size();
	for (int Set = 0; Set < NumberOfSets; Set++) {
		for (int index = 0; index < S.size(); ++index) {
			//Return_Vector.at(Set*S.size() + index) = S.at(index) + Sc.at(Set);
			Return_Vector.push_back(S.at(index) + Sc.at(Set));
		}
		Return_Vector.push_back(CHAR_MAX);
	}
	return Return_Vector;
}

void LDPC_Code_16_12(MATRIX<__int8> &H) {
	vector<__int8> Z = { 0,2,1,3 };
	for (int h_matrix = 0; h_matrix < 4; h_matrix++) {
		for (int j = 0; j < 4; j++) {
			int col = (Z.at(h_matrix) + j) % 4;
			H._matrix[j][col + h_matrix * 4] = 1;
		}
	}
}

void Convert_SystematicG_to_H(MATRIX<__int8> &H, MATRIX<__int8> &G){
	int parity_num = G.Col_number - G.Row_number;
	H.Building_Empty_Matrix(parity_num, G.Col_number);
	for (int row = 0; row < parity_num; row++) {
		for (int col = 0; col < G.Row_number; col++) {
			H._matrix[row][col] = G._matrix[col][row + G.Row_number];
		}
	}
	for (int pivot = 0; pivot < parity_num; pivot++) {
		H._matrix[pivot][pivot + G.Row_number] = 1;
	}
}

inline void Matrix_Mul(MATRIX<__int8> & M1, MATRIX<__int8> & M2, MATRIX<__int8>& res) {
	if (M1.Col_number != M2.Row_number || !(M1.Row_number == res.Row_number)
		|| !(M2.Col_number == res.Col_number)) {
		cout << "Matrix_Mul error";
		return;
	}
	int tmp;
	for (int i = 0; i < M1.Row_number; i++) {
		for (int j = 0; j < M2.Col_number; j++) {
			tmp = 0;
			for (int k = 0; k < M1.Col_number; k++) {
				tmp ^= M1._matrix[i][k] * M2._matrix[k][j];
			}
			res._matrix[i][j] = tmp;
		}
	}
	
}
inline void GH_test(MATRIX<__int8> & M1, MATRIX<__int8> & M2) {
	if (M1.Col_number != M2.Col_number ) {
		cout << "Matrix_test error";
		return;
	}
	int tmp;
	for (int i = 0; i < M1.Row_number; i++) {
		for (int j = 0; j < M2.Row_number; j++) {
			tmp = 0;
			for (int k = 0; k < M1.Col_number; k++) {
				tmp ^= M1._matrix[i][k] * M2._matrix[j][k];
			}
			if (tmp != 0) {
				cout << "G H isn't pair" << endl;
				return;
			}
		}
	}
	cout << "correct G H pair" << endl;

}

inline void repeatition(MATRIX<__int8> & G) {
	G.Col_number *= 2;
	for (int i = 0; i < G.Row_number; i++) {
		G._matrix.at(i).resize(G.Col_number);
	}
	for (int i = 0; i < G.Row_number; i++) {
		for (int j = 0; j < G.Col_number / 2; j++) {
			G._matrix.at(i).at(j + G.Col_number / 2) = G._matrix.at(i).at(j);
		}
	}
}

inline void convert_H_G(MATRIX<__int8> & H,MATRIX<__int8> & G) {
	//G identity sub matrix
	for (int i = 0; i < G.Row_number; i++) {
		G._matrix[i][i] = 1;
	}
	//transform P -> P^T
	for (int i = 0; i < H.Row_number; i++) {
		for (int j = 0; j < G.Row_number; j++) {
			G._matrix[j][G.Row_number + i] = H._matrix[i][j];
		}
	}
}

inline void exchang_column(MATRIX<__int8> & G, int level) {
	MATRIX<__int8> tmp = G;
	for (int i = level; i < G.Col_number; i++) {
		for (int j = 0; j < G.Row_number; j++) {
			G._matrix[j][i - level] = tmp._matrix[j][i];
		}
	}
	int start = G.Col_number - level;
	for (int i = 0; i < level; i++) {
		for (int j = 0; j < G.Row_number; j++) {
			G._matrix[j][start + i] = tmp._matrix[j][i];
		}
	}
}
