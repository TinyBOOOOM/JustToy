#include "RTrees.h"

#define HAND_INPUT

int main( int argc, char** argv ) {
#ifdef HAND_INPUT
    int year;
    std::string eddyfile;
    std::string argofile;
    std::string matchfile;

    std::cout << "Year: ";
    std::cin >> year;
    std::cout << "Eddy file path: ";
    std::cin >> eddyfile;
    std::cout << "Argo file file path: ";
    std::cin >> argofile;
    std::cout << "Match result file path: ";
    std::cin >> matchfile;

    RTreeMatchTest( year, eddyfile, argofile, matchfile );
#else
    RTreeMatchTest( 2001, "D:/Projects/Quad-System/Quad-System/eddy_global_2001", "argo_all_2001.txt", "Match_RES.txt" );
#endif
	system("pause");
    return 0;
}
