//Functions for flooding digital elevation models
//Based on source code in RichDEM (Barnes, Richard. 2016. RichDEM: Terrain Analysis Software. http://github.com/r-barnes/richdem)
//Modified to Rcpp by Kenneth Thorø Martinsen

#include <Rcpp.h>
//#include <vector>
#include <queue>
//#include <cmath>
//#include <functional>
using namespace Rcpp;
using namespace std;

class cell{
  public:
    int r;
    int c;
    cell(){}
    cell(int r, int c): r(r), c(c){}
};

class cellz: 
  public cell{
  public:
    double z;
    cellz(){}
    cellz(int r, int c, double z): cell(r, c), z(z){}
    bool operator> (const cellz &a) const {return z > a.z; }
};

// Flow coordinate system
// 234
// 105
// 876

//Facet                0   1   2   3   4  5  6  7   8
const int d8[9] =     {0,  1,  2,  3,  4, 5, 6, 7,  8};
const int d8_inv[9] = {0,  5,  6,  7,  8, 1, 2, 3,  4};
const int dx[9] =     {0, -1, -1,  0,  1, 1, 1, 0, -1}; ///< x offsets of D8 neighbours, from a central cell
const int dy[9] =     {0,  0, -1, -1, -1, 0, 1, 1,  1}; ///< y offsets of D8 neighbours, from a central cell

// Constants used in functions
int labels_nodata = 0;
double dem_nodata = -9999.0;
double pittop = dem_nodata;
int false_pit_cells = 0;
int clabel = 1;

// [[Rcpp::export]]
NumericMatrix pf_barnes2014(NumericMatrix dem){
  
  // Improved priority flood (algorithm 2) in:
  // "Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024"
  
  // Initialize structures and counters
  priority_queue<cellz, vector<cellz>, greater<cellz>> open;
  queue<cellz> pit;
  LogicalMatrix closed(dem.nrow(), dem.ncol());
  
  // Add edge cells to priority queue
  // Horizontal edges
  for(int x = 0; x < dem.ncol(); x++){
    open.push(cellz(0, x, dem(0, x)));
    open.push(cellz(dem.nrow()-1, x, dem(dem.nrow()-1, x)));
    closed(0, x) = true;
    closed(dem.nrow()-1, x) = true;
  }
  
  // Vertical edges
  for(int y = 0; y < dem.nrow(); y++){
    open.push(cellz(y, 0, dem(y, 0)));
    open.push(cellz(y, dem.ncol()-1, dem(y, dem.ncol()-1)));
    closed(y, 0) = true;
    closed(y, dem.ncol()-1) = true;
  }
  
  while(open.size()>0 || pit.size()>0){
    cellz c;
    if(pit.size()>0){
      c = pit.front();
      pit.pop();
    } else {
      c=open.top();
      open.pop();
    }

    for(int n=1; n<=8; n++){
      int nc=c.c+dx[n];
      int nr=c.r+dy[n];
      if(!(0<=nc && nc<dem.ncol() && 0<=nr && nr<dem.nrow())) 
        continue;
      if(closed(nr,nc))
        continue;
      
      closed(nr,nc) = true;
      
      if(dem(nr,nc) <= c.z){
        
        if(dem(nr,nc) < c.z){
          
          dem(nr,nc) = c.z;
          
        }
        
        pit.push(cellz(nr,nc,c.z));
        
      } else {
        open.push(cellz(nr, nc, dem(nr, nc)));
        }
      }
    }
  
  return(dem);
  
}

// [[Rcpp::export]]
NumericMatrix pf_eps_barnes2014(NumericMatrix dem){
  
  // Improved priority flood (algorithm 3) in:
  // "Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024"
  
  // Initialize structures and counters
  priority_queue<cellz, vector<cellz>, greater<cellz>> open;
  queue<cellz> pit;
  LogicalMatrix closed(dem.nrow(), dem.ncol());
  
  // Add edge cells to priority queue
  // Horizontal edges
  for(int x = 0; x < dem.ncol(); x++){
    open.push(cellz(0, x, dem(0, x)));
    open.push(cellz(dem.nrow()-1, x, dem(dem.nrow()-1, x)));
    closed(0, x) = true;
    closed(dem.nrow()-1, x) = true;
  }
  
  // Vertical edges
  for(int y = 0; y < dem.nrow(); y++){
    open.push(cellz(y, 0, dem(y, 0)));
    open.push(cellz(y, dem.ncol()-1, dem(y, dem.ncol()-1)));
    closed(y, 0) = true;
    closed(y, dem.ncol()-1) = true;
  }
  
  while(open.size()>0 || pit.size()>0){
    cellz c;
    
    if(pit.size()>0 && open.size()>0 && open.top().z == pit.front().z){
      c = open.top();
      open.pop();
      pittop = dem_nodata;
    } else if(pit.size()>0){
      c=pit.front();
      pit.pop();
      if(pittop == dem_nodata)
        pittop = dem(c.r,c.c);
    } else {
      c=open.top();
      open.pop();
      pittop= dem_nodata;
    }
    
    for(int n=1; n<=8; n++){
      
      int nc=c.c+dx[n];
      int nr=c.r+dy[n];
      
      if(!(0<=nc && nc<dem.ncol() && 0<=nr && nr<dem.nrow())) 
        continue;
      
      if(closed(nr,nc))
        continue;
      
      closed(nr,nc) = true;
      
      if(dem(nr,nc) == dem_nodata){
        pit.push(cellz(nr, nc, dem_nodata));
      }
      
      else if(dem(nr, nc) <= nextafter(c.z, numeric_limits<double>::infinity())){
        if(pittop != dem_nodata && pittop < dem(nr,nc) && nextafter(c.z, numeric_limits<double>::infinity()) >= dem(nr,nc))
          ++false_pit_cells;
        dem(nr,nc) = nextafter(c.z, numeric_limits<double>::infinity());
        pit.push(cellz(nr, nc, dem(nr, nc)));
      } else
        open.push(cellz(nr, nc, dem(nr, nc)));
    }
  }
  
  if(false_pit_cells)
    Rcout<<"Warning: While raising elevation of depression cells. Elevation for " <<false_pit_cells<< "  cells were increased above that of surrounding cells." << std::endl;
  
  return(dem);
}


// [[Rcpp::export]]
List pf_watersheds_barnes2014(NumericMatrix dem){
  
  // Improved priority flood with watershed labels (algorithm 5) in:
  // "Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024"
  
  // Initialize structures and counters
  priority_queue<cellz, vector<cellz>, greater<cellz>> open;
  queue<cellz> pit;
  LogicalMatrix closed(dem.nrow(), dem.ncol());
  IntegerMatrix labels(dem.nrow(), dem.ncol());
  
  // Add edge cells to priority queue
  // Horizontal edges
  for(int x = 0; x < dem.ncol(); x++){
    open.push(cellz(0, x, dem(0, x)));
    open.push(cellz(dem.nrow()-1, x, dem(dem.nrow()-1, x)));
    closed(0, x) = true;
    closed(dem.nrow()-1, x) = true;
  }
  
  // Vertical edges
  for(int y = 0; y < dem.nrow(); y++){
    open.push(cellz(y, 0, dem(y, 0)));
    open.push(cellz(y, dem.ncol()-1, dem(y, dem.ncol()-1)));
    closed(y, 0) = true;
    closed(y, dem.ncol()-1) = true;
  }
  
  while(open.size()>0 || pit.size()>0){
    cellz c;
    if(pit.size()>0){
      c = pit.front();
      pit.pop();
    } else {
      c=open.top();
      open.pop();
    }
    
    if(labels(c.r,c.c) == labels_nodata && dem(c.r,c.c) != dem_nodata){
      labels(c.r,c.c) = clabel++;
      }
      
    for(int n=1; n<=8; n++){
      int nc=c.c+dx[n];
      int nr=c.r+dy[n];
      if(!(0<=nc && nc<dem.ncol() && 0<=nr && nr<dem.nrow())) 
        continue;
      if(closed(nr,nc))
        continue;
      
      labels(nr,nc) = labels(c.r,c.c);
      
      closed(nr,nc) = true;
      
      if(dem(nr,nc) <= c.z){
        
        if(dem(nr,nc) < c.z){
          
          dem(nr,nc) = c.z;
          
        }
        
        pit.push(cellz(nr,nc,c.z));
        
      } else {
        
        open.push(cellz(nr, nc, dem(nr, nc)));
        
      }
    }
  }
  
  List result = List::create(_["dem"] = dem, _["labels"] = labels);
  
  return(result);
  
}
