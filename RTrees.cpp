#include "RTrees.h"

RTreeAnnual::RTreeAnnual( const std::string & path, int year ) : isInit_( readEddyFile( path, year ) ) {}

RTreeAnnual::RTreeAnnual() : isInit_( false ) {}

bool RTreeAnnual::readEddyFile( const std::string & path, int year ) {
    eddys_.clear();
    eddys_.reserve( EDDY_PER_YEAR_MAX );

    std::cout << "Open path: " << path << std::endl;

    vector<char*> buffers( DAYS_PER_YEAR_MAX, NULL );
    {
        ShowTime st( "Read eddy file" );
        unsigned int dayBeg = year * 1000 + 1;
        for ( size_t i = 0; i < DAYS_PER_YEAR_MAX; ++i ) {
            std::string filename( to_string( dayBeg++ ) + ".txt" );
            ifstream ifs( path + "/" + filename, ios::binary );
            if ( ifs.is_open() ) {
                std::cout << "\tRead file: " << i + 1 << "/" << DAYS_PER_YEAR_MAX << "\r";
                // 一次性读入内存
                filebuf * pbuffer = ifs.rdbuf();
                unsigned int size = ( unsigned int )pbuffer->pubseekoff( 0, ios::end, ios::in );
                pbuffer->pubseekpos( 0, ios::in );
                char * buffer = new char[size];
                pbuffer->sgetn( buffer, size );
                buffers[i] = buffer;

                ifs.close();
            }
        }
        std::cout << "\n";
    }

    // 解析eddy数据
    {
        ShowTime st( "Parse eddy data" );
        unsigned int eddyNum = 0;
        for ( size_t i = 0; i < DAYS_PER_YEAR_MAX; ++i ) { // 365天的数据
            if ( NULL != buffers[i] ) {
                std::cout << "\tDeal process: " << i + 1 << "/" << DAYS_PER_YEAR_MAX << "\r";

                stringstream sss( buffers[i] );
                char buffff[64] = { 0 }; // 行缓冲区
                while ( sss.getline( buffff, 64 ) ) {
                    Eddy_s eddy = { 0 };
                    // 格式化输入
                    sscanf_s( buffff, "%f, %f, %f, %f, %f,"
                              , &eddy.lon, &eddy.lat, &eddy.p1, &eddy.p2, &eddy.p3 );
                    eddy.id = eddyNum; // 自增序号
                    eddy.day = i + 1; // 记录天数
                    eddys_.push_back( eddy );
                    ++eddyNum;
                }
            }
        }
        std::cout << "\nTotal eddy number: " << eddyNum << std::endl;
    }

    // 清空buffer
    for ( auto iter = buffers.begin(); iter != buffers.end(); ++iter ) {
        delete[] * iter;
        *iter = NULL;
    }
    return true;
}

void RTreeAnnual::creatTree() {
    if ( eddys_.empty() ) {
        isInit_ = false;
        return;
    }
    rtrees_.clear();
    rtrees_.resize( DAYS_PER_YEAR_MAX );
    ShowTime st( std::string( "Creat R-Tree with " + to_string( eddys_.size() ) + " eddy data" ) );
    for ( size_t i = 0; i < eddys_.size(); ++i ) {
        Eddy_s & eddy = eddys_[i];
        box_t b( point_t( eddy.lon - BOX_RADIUS, eddy.lat - BOX_RADIUS ), point_t( eddy.lon + BOX_RADIUS, eddy.lat + BOX_RADIUS ) );
        rtrees_[eddy.day - 1].insert( make_pair( b, eddy.id ) );
    }
    isInit_ = true;
}

bool RTreeAnnual::query( int day, float lon, float lat, float radius, std::vector<value_t> & results ) {
    if ( !isInit_ )
        return false;
    if ( !( 0 < day&&day <= 366 ) )
        return false;

    box_t b( point_t( lon - radius, lat - radius ), point_t( lon + radius, lat + radius ) );
    rtrees_[day - 1].query( bgi::intersects( b ), std::back_inserter( results ) );
    return true;
}

bool readArgoFile( const std::string & filepath, vector<Argo_s> & argos, unsigned int max_data /*= 100000U*/ ) {
    ShowTime st( "Read argo file" ); // 读取并处理argo数据
    ifstream ifs( filepath );
    if ( !ifs.is_open() )
        return false;
    argos.clear();
    argos.reserve( max_data ); // 分配设置好的内存量
    while ( !ifs.eof() ) {
        char bufffer[64] = { 0 };
        ifs.getline( bufffer, 64 );
        Argo_s argo = { 0 };
        float idTmp = 0.0f; // 输入指数形式的数据，以float形式解析，而后转uint
        sscanf_s( bufffer, "%f, %f, %d, %d, %f,"
                  , &argo.lon, &argo.lat, &argo.day, &argo.radius, &idTmp );
        argo.id = ( unsigned int )idTmp;
        argos.push_back( argo );
    }
    return true;
}

bool writeMatchFile( const vector<std::pair<unsigned int, unsigned int>> & matchPairs, const RTreeAnnual & rt, const vector<Argo_s> & argos, const std::string & filepath ) {
    ShowTime outst( "Write result" );
    if ( matchPairs.empty() )
        return false;

    ofstream ofs( filepath, ios::trunc );
    if ( !ofs.is_open() )
        return false;
    // 匹配数量
    ofs << "Match number: " << matchPairs.size() << std::endl;
    // 表头
    ofs << "eddy_lon eddy_lat eddy_radius eddy_zhenfu eddy_v argo_lon argo_lat argo_day argo_radius argo_id\n";
    for ( size_t i = 0; i < matchPairs.size(); ++i ) {
        const Eddy_s & eddy = rt.at( matchPairs[i].first );
        const Argo_s & argo = argos[matchPairs[i].second];
        ofs << eddy.lon << ' ' << eddy.lat << ' ' << eddy.p1 << ' ' << eddy.p2 << ' ' << eddy.p3 << ' ';
        ofs << argo.lon << ' ' << argo.lat << ' ' << argo.day << ' ' << argo.radius << ' ' << argo.id << '\n';
    }
    ofs.close();
    return true;
}

void RTreeMatchTest( int year, const std::string & eddyfile, const std::string & argofile, const std::string & outfile ) {
    ShowTime allst( "ALL TIME" );

    // 建树
    RTreeAnnual rt;
    // 先读数据
    rt.readEddyFile( eddyfile, year );
    // 再建R-Tree
    rt.creatTree();

    MatchPairs_t matchPairs; // 匹配数据
    matchPairs.reserve( 500000 ); // 预先设置50W数据的内存

    // 读取并处理argo数据
    vector<Argo_s> argos;
    readArgoFile( argofile, argos, 100000U );
    std::cout << "Argo number: " << argos.size() << std::endl;

    // 匹配
    {
        ShowTime st( "Match eddy and argo" );
        unsigned long matchNum = 0;

        for ( size_t i = 0; i < argos.size(); ++i ) {
            vector<RTreeAnnual::value_t> res; // 存单次匹配结果
            rt.query( argos[i].day, argos[i].lon, argos[i].lat, 1.5f, res ); // 取1.5的半径
            if ( !res.empty() ) {
                matchNum += res.size();
                for ( size_t j = 0; j < res.size(); ++j ) {
                    // 插入每对匹配数据
                    matchPairs.push_back( std::make_pair( res[j].second, i ) ); // <eddy_id, argo_id>
                }
            }
        }

        std::cout << "Match number: " << matchNum << std::endl;
    }

    // 写入文件
    writeMatchFile( matchPairs, rt, argos, outfile );
}
