/************************************************************************/
/* 采用R-Tree作为基本空间索引												*/
/* 每年中的每一天都创建一颗R-Tree											*/
/* 每个点看成一个正方形，可设置其对角线半径									*/
/* 以eddy数据为数据集，argo数据作匹配										*/
/************************************************************************/

#ifndef _RTREES_H_
#define _RTREES_H_

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>

// to store queries results
#include <vector>

// i/o
#include <iostream>
#include <fstream>

#include <time.h>

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// 经度，纬度，半径，振幅，最大地转流速度，某天，自定义编号
struct Eddy_s {
    float lon;
    float lat;
    float p1;
    float p2;
    float p3;
    int day;
    unsigned int id;
};

struct Argo_s {
    float lon;
    float lat;
    int day;
    int radius;
    unsigned int id;
};

class ShowTime { // 用生命周期控制显示时间
  public:
    explicit ShowTime( const std::string & name )
        : time_( clock() ), name_( name ) {
		std::cout << "----------------------------------------\n";
		std::cout << name_ << " begin\n";
		
    }
    ~ShowTime() {
        std::cout << name_ << " spend time: " << clock() - time_ << " ms\n";
		std::cout << "----------------------------------------\n";
    }
  private:
    clock_t time_;
    std::string name_;
};

class RTreeAnnual {  // 365(366)天，每天一个R-Tree
  public:
    const unsigned int DAYS_PER_YEAR_MAX = 366;
    const float BOX_RADIUS = 1.5f;
    const unsigned int DATA_PER_DAY_MAX = 5000;
    const unsigned int EDDY_PER_YEAR_MAX = 10000000;

    typedef bg::model::point<float, 2, bg::cs::cartesian> point_t;
    typedef bg::model::box<point_t> box_t;
    typedef std::pair<box_t, unsigned> value_t;
    typedef bgi::rtree<value_t, bgi::quadratic<16>> rtree_t;
    typedef std::vector<rtree_t> rtrees_t;
    typedef std::vector<Eddy_s> eddys_t;

  public:
    RTreeAnnual();
    RTreeAnnual( const std::string & path, int year );
    ~RTreeAnnual() {};
    // 读取并解析eddy文件
    bool readEddyFile( const std::string & path, int year );
    // 创建树节点
    void creatTree();
    // radius单位是经纬度，几何意义是正方形的对角线半径
    bool query( int day, float lon, float lat, float radius, std::vector<value_t> & results );
    inline Eddy_s at( unsigned int index ) const {
        return eddys_[index];
    };
    inline Eddy_s at( unsigned int index ) {
        return eddys_[index];
    };

  private:
    bool isInit_;
    rtrees_t rtrees_;
    eddys_t eddys_;
};

// filepath - 路径/文件
// argos - 接收容器
// max_data - 最大数据量（只要比总数据量大，就能加快速度，默认10万）
bool readArgoFile( const std::string & filepath, vector<Argo_s> & argos, unsigned int max_data = 100000U );

typedef vector<std::pair<unsigned int, unsigned int>> MatchPairs_t;
// 写入文件
bool writeMatchFile( const MatchPairs_t & matchPairs, const RTreeAnnual & rt, const vector<Argo_s> & argos, const std::string & filepath );
// 测试匹配
void RTreeMatchTest( int year, const std::string & eddyfile, const std::string & argofile, const std::string & outfile );

#endif // _RTREES_H_
