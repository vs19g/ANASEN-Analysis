#ifndef Mapping_h
#define Mapping_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <utility>

#include <TMath.h>

const std::map<int, unsigned short> board = {
   {0, 17122},  // id, sn
   {1, 17123},
   {2, 22320},
   {3, 22130},
   {4, 22129},
   {5, 15529},
   {6, 15528},
   {7,   379},
   {8,   409},
   {9,   405}
};
const int nBd = board.size();

const int nV1740 = 7;
const int nV1725 = 3;

//+++++++++++++++++++ detID;
// The detectors are seperated into 2 type: SuperX3, QQQ, and PC
// the SuperX3 has 24 detectors for each kind, wach detector has 12 channels
// the QQQ has 4 detectors for each kind, each detector has 32 channels
// the PC has 2 types, anode and cathode, each has 24 channels
// The detID = Type * 10000 + index * 100 + channel
// fro example, detID(superX3-8, ch-7) = 00807


// use the GenMapping() to get that 
const std::vector<int> mapping = {
  
  //================== 17122
   206,   207,   204,   205,   203,   202,   201,   200,   406,   407,   404,   405,   403,   402,   401,   400,
     6,     7,     4,     5,     3,     2,     1,     0,   506,   507,   504,   505,   503,   502,   501,   500,
   111,   110,   109,   108,   211,   210,   209,   208,   311,   310,   309,   308,   411,   410,   409,   408,
   106,   107,   104,   105,   103,   102,   101,   100,   306,   307,   304,   305,   303,   302,   301,   300,
  //================== 17123
   606,   607,   604,   605,   603,   602,   601,   600,  1106,  1107,  1104,  1105,  1103,  1102,  1101,  1100,
   711,   710,   709,   708,   811,   810,   809,   808,   911,   910,   909,   908,  1011,  1010,  1009,  1008,
   706,   707,   704,   705,   703,   702,   701,   700,   906,   907,   904,   905,   903,   902,   901,   900,
   806,   807,   804,   805,   803,   802,   801,   800,  1006,  1007,  1004,  1005,  1003,  1002,  1001,  1000,
  //================== 22320
  1911,  1910,  1909,  1908,  2011,  2010,  2009,  2008,  2111,  2110,  2109,  2108,  2211,  2210,  2209,  2208,
  1906,  1907,  1904,  1905,  1903,  1902,  1901,  1900,  2106,  2107,  2104,  2105,  2103,  2102,  2101,  2100,
  1806,  1807,  1804,  1805,  1803,  1802,  1801,  1800,  2306,  2307,  2304,  2305,  2303,  2302,  2301,  2300,
  2006,  2007,  2004,  2005,  2003,  2002,  2001,  2000,  2206,  2207,  2204,  2205,  2203,  2202,  2201,  2200,
  //================== 22130
  1311,  1310,  1309,  1308,  1411,  1410,  1409,  1408,  1511,  1510,  1509,  1508,  1611,  1610,  1609,  1608,
    11,    10,     9,     8,   511,   510,   509,   508,   611,   610,   609,   608,  1111,  1110,  1109,  1108,
  1406,  1407,  1404,  1405,  1403,  1402,  1401,  1400,  1606,  1607,  1604,  1605,  1603,  1602,  1601,  1600,
  1306,  1307,  1304,  1305,  1303,  1302,  1301,  1300,  1506,  1507,  1504,  1505,  1503,  1502,  1501,  1500,
  //================== 22129
 10015, 10014, 10013, 10012, 10011, 10010, 10009, 10008, 10007, 10006, 10005, 10004, 10003, 10002, 10001, 10000,
  1206,  1207,  1204,  1205,  1203,  1202,  1201,  1200,  1706,  1707,  1704,  1705,  1703,  1702,  1701,  1700,
 10115, 10114, 10113, 10112, 10111, 10110, 10109, 10108, 10107, 10106, 10105, 10104, 10103, 10102, 10101, 10100,
 10116, 10117, 10118, 10119, 10120, 10121, 10122, 10123, 10124, 10125, 10126, 10127, 10128, 10129, 10130, 10131,
  //================== 15529
 10016, 10017, 10018, 10019, 10020, 10021, 10022, 10023, 10024, 10025, 10026, 10027, 10028, 10029, 10030, 10031,
 10215, 10214, 10213, 10212, 10211, 10210, 10209, 10208, 10207, 10206, 10205, 10204, 10203, 10202, 10201, 10200,
 10216, 10217, 10218, 10219, 10220, 10221, 10222, 10223, 10224, 10225, 10226, 10227, 10228, 10229, 10230, 10231,
 10315, 10314, 10313, 10312, 10311, 10310, 10309, 10308, 10307, 10306, 10305, 10304, 10303, 10302, 10301, 10300,
  //================== 15528
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
 10316, 10317, 10318, 10319, 10320, 10321, 10322, 10323, 10324, 10325, 10326, 10327, 10328, 10329, 10330, 10331,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
  //================== 379
 20116, 20117, 20118, 20119, 20120, 20121, 20122, 20123, 20016, 20017, 20018, 20019, 20020, 20021, 20022, 20023,
  //================== 409
 20000, 20001, 20002, 20003, 20004, 20005, 20006, 20007, 20008, 20009, 20010, 20011, 20012, 20013, 20014, 20015,
  //================== 405
 20100, 20101, 20102, 20103, 20104, 20105, 20106, 20107, 20108, 20109, 20110, 20111, 20112, 20113, 20114, 20115,
  
};

void PrintMapping(){

  int digiID = 0;
  int count = 0;
  printf("==================== ID-MAP: \n");
  printf("%11s|", "");  for(int i = 0 ; i < 16; i++ ) printf("%7d|", i);
  printf("\n");
  for(int i = 0 ; i < 12 + 16*8; i++ ) printf("-");
  for(size_t i = 0 ; i < mapping.size(); i ++){
    if( (i) % 16 == 0 ) {
      printf("\n");

      if( digiID < nBd ){
        if( board.at(digiID) > 1000 ) {
          if( count == 3 ) digiID ++;
          if( i % 64 == 0 ) {
            printf("%11d|", board.at(digiID));
            count = 0;
          }
        }else{
          if( count == 1 ) digiID ++;
          if( i % 16 == 0 ) {
            printf("%11d|", board.at(digiID));
            count = 0;
          }
        }
      }

      if( count != 0) printf("%11s|", "");
      count ++;
    }

    int typeID = mapping[i] / 10000;
    int detID = (mapping[i] - typeID*10000 )/100;
    int ch = mapping[i] - typeID*10000 - detID * 100;

    if( mapping[i] == -1 ) {

      printf("%7s|", "");

    }else{

      if( typeID == 0){ // SX3

        printf("\033[36m%3d(%2d)\033[0m|", detID, ch); 

      }else if( typeID == 1){ // QQQ

        printf("\033[91m%3d(%2d)\033[0m|", detID, ch);
        
      }else if( typeID == 2){ // PC

        printf("\033[35m%3d(%2d)\033[0m|", detID, ch);
        
      }else{


      }
    }
  }
  printf("\n");
  for(int i = 0 ; i < 12 + 16*8; i++ ) printf("-");
  printf("\n");

}


void GenMapping(std::string mapFile){


  std::vector<int> map;


  std::ifstream inputFile(mapFile);  // Replace "your_file.txt" with the actual file path

  if (!inputFile.is_open()) {
    printf("Error: Could not open the file (%s).\n", mapFile.c_str());
    return ;
  }

  std::string line;

  // Read the file line by line
  while (std::getline(inputFile, line)) {
    std::vector<std::string> words;
    std::istringstream iss(line);

    // Extract words from the current line
    while (true) {
      std::string word;
      if (!(iss >> word)) break;  // Break if there are no more words

      word.erase(std::remove_if(word.begin(), word.end(), ::isspace), word.end());
      words.push_back(word);
    
    }

    if( atoi(words[0].c_str()) % 16 == 0 ) printf("=================\n");

    for( size_t i = 0; i < words.size(); i++) printf("|%9s", words[i].c_str());

    int detID = atoi(words[1].c_str())*100;

    if( words[2] == "BARREL" ) {
      if( words[3] == "FRONTDOWN" ){
        int chID = atoi(words[4].c_str());
        if( chID % 2 != 0 ) chID -= 1;
        detID += chID;
      }

      if( words[3] == "FRONTUP" ){
        int chID = atoi(words[4].c_str());
        if( chID % 2 == 0 ) chID += 1;
        detID += chID;
      }

      if( words[3] == "BACK") detID += atoi(words[4].c_str()) + 8;
    }

    if( words[2] == "FQQQ" ) {
      detID += 10000;
      if( words[3] == "WEDGE") detID += atoi(words[4].c_str());
      if( words[3] == "RING") detID += atoi(words[4].c_str()) + 16;
    }

    if( words[2] == "PC" ) {
      detID += 20000;
      if( words[3] == "ANODE") detID += atoi(words[4].c_str());
      if( words[3] == "CATHODE") detID += 100 + atoi(words[4].c_str());

    }

    if( words[2] == "blank") {
      detID = -1;
    }

    map.push_back(detID);

    printf("|  %5d", detID);
    printf("|\n");
  }

  // Close the file
  inputFile.close();

  int digiID = 0;
  int count = 0;
  printf("===============================\n");
  for( size_t i = 0; i < ((map.size() +15)/16) * 16; i++){
    if( i % 16 == 0) {
      printf("\n");
      if( digiID < nBd ){
        if( board.at(digiID) > 1000 ) {
          if( count == 3 ) digiID ++;
          if( i % 64 == 0 ) {
            printf("  //================== %d\n", board.at(digiID));
            count = 0;
          }
        }else{
          if( count == 1 ) digiID ++;
          if( i % 16 == 0 ) {
            printf("  //================== %d\n", board.at(digiID));
            count = 0;
          }
        }
      }
      count ++;
    }
    if( i < map.size() ){
      printf(" %5d,", map[i]);
    }else{
      printf(" %5d,", -1);
    }
  }
  printf("\n\n===============================\n");

  printf("sorting mapping and see if there any repeated\n");
  std::sort(map.begin(), map.end());

  for( size_t i = 1; i < map.size(); i++){
    if( map[i] == -1 ) continue;
    if( map[i] == map[i-1] ) printf("%5d \n", map[i]);
  }
  printf("=========== Done. if nothing show, no repeat. \n");


}

#endif 