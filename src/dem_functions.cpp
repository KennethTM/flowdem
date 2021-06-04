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
int watershed_nodata = 0;
int flowdir_nodata = 0;
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
  //bool closed[dem.nrow()][dem.ncol()];
  
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
List pf_basins_barnes2014(NumericMatrix dem){
  
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

static int d8_flowdir(NumericMatrix &dem, const int r, const int c){
  
  // Helper function for determining d8 flow directions (RichDEM)
  
  double minimum_elevation = dem(r, c);
  int flowdir = 0;
  
  // Flow direction for edge cells
  if (r==0 || c==0 || c == dem.ncol()-1 || r == dem.nrow()-1){
    if(r==0 && c==0)
      return 2;
    else if(r == dem.nrow()-1 && c == 0)
      return 8;
    else if(r == 0 && c == dem.ncol()-1)
      return 4;
    else if(r == dem.nrow()-1 && c == dem.ncol()-1)
      return 6;
    else if(c == 0)
      return 1;
    else if(c == dem.ncol()-1)
      return 5;
    else if(r == 0)
      return 3;
    else if(r == dem.nrow()-1)
      return 7;
  }
  
  // Flow direction for all other cells
  for(int n=1; n<=8; n++){
    if(dem(r+dy[n], c+dx[n]) < minimum_elevation || (dem(r+dy[n], c+dx[n]) == minimum_elevation && flowdir>0 && flowdir%2==0 && n%2==1)){
      minimum_elevation = dem(r+dy[n], c+dx[n]);
      flowdir = n;
    }
  }
  
  return flowdir;
  
}

// [[Rcpp::export]]
IntegerMatrix d8_flow_directions(NumericMatrix dem){
  
  // Function for determining d8 flow directions (RichDEM)
  
  IntegerMatrix flowdirs(dem.nrow(), dem.ncol());

  //#pragma omp parallel for
  //see https://mfasiolo.github.io/sc2-2019/rcpp_advanced_iii/2_rcppparallel/
  for(int r=0; r < dem.nrow(); r++){
    for(int c=0; c < dem.ncol(); c++)
      if(dem(r, c) == dem_nodata)
        flowdirs(r, c) = flowdir_nodata;
      else
        flowdirs(r, c) = d8_flowdir(dem, r, c);
  }
  
  return flowdirs;
}


// [[Rcpp::export]]
NumericMatrix d8_flow_accum(IntegerMatrix flowdirs){
  
  // Function for determining d8 flow accumulation (RichDEM)
  
  std::queue<cell> sources;
  
  NumericMatrix dependency(flowdirs.nrow(), flowdirs.ncol());
  NumericMatrix area(flowdirs.nrow(), flowdirs.ncol());
  double area_nodata = -1;
  
  //#pragma omp parallel for
  for(int r = 0; r<flowdirs.nrow(); r++){
    for(int c = 0; c<flowdirs.ncol(); c++){
      if(flowdirs(r, c) == flowdir_nodata){
        area(r,c) = area_nodata;
        continue;
      }
      
      int n = flowdirs(r,c); 
      if(n == flowdir_nodata) //check if correct (NOFLOW)
        continue;
      
      int nc = c+dx[n];
      int nr = r+dy[n];
      
      if(!(0<=nc && nc<flowdirs.ncol() && 0<=nr && nr<flowdirs.nrow()))
        continue;
      
      ++dependency(nr,nc);
    }
  }

  for(int r=0; r<flowdirs.nrow(); r++){
    for(int c = 0; c<flowdirs.ncol(); c++){
      if(dependency(r,c) == 0 && flowdirs(r, c) != flowdir_nodata)
        sources.emplace(cell(r, c));
    }
  }
      
  long int ccount=0;
  while(sources.size()>0){
    cell c = sources.front();
    sources.pop();
    
    ccount++;

    area(c.r,c.c)++;
    
    int n = flowdirs(c.r,c.c);
    
    if(n == flowdir_nodata)
      continue;
    
    int nc = c.c+dx[n];
    int nr = c.r+dy[n];
    
    if(!(0<=nc && nc<flowdirs.ncol() && 0<=nr && nr<flowdirs.nrow()))
      continue;
    
    if(flowdirs(nr, nc) == flowdir_nodata)
      continue;
    
    area(nr,nc) += area(c.r,c.c);
    --dependency(nr,nc);
        
    if(dependency(nr,nc) == 0)
      sources.emplace(cell(nr, nc));
  }

  return area;

}

// // [[Rcpp::export]]
// IntegerMatrix d8_watershed(IntegerMatrix flowdirs, NumericMatrix target_rc){
//   
//   // Function for d8 watersheds to a target area identified by row-col indexes
//   IntegerMatrix watershed(flowdirs.nrow(), flowdirs.ncol());
//   
//   std::queue<cell> expansion;
// 
//   for(int r = 0; r < target_rc.nrow(); r++){
//     expansion.push(cell(target_rc(r, 0)-1, target_rc(r, 1)-1));
//   }
//   
//   while(expansion.size()>0){
//     cell c = expansion.front();
//     expansion.pop();
//     
//     for(int n = 1; n <= 8; n++){
//       
//       int nc = c.c+dx[n];
//       int nr = c.r+dy[n];
//       
//       if(!(0<=nc && nc<flowdirs.ncol() && 0<=nr && nr<flowdirs.nrow()))
//         continue;
//       
//       else if(flowdirs(nr, nc) == flowdir_nodata)
//         continue;
//       
//       else if(watershed(nr, nc) == watershed_nodata && n == d8_inv[flowdirs(nr, nc)]){
//         expansion.push(cell(nr, nc));
//         watershed(nr, nc) = 1;
//       }
//     }
//   }
//   
//   return watershed;
// }

// [[Rcpp::export]]
IntegerMatrix d8_watershed_nested(IntegerMatrix flowdirs, NumericMatrix target_rc, bool nested){
  
  // Function for d8 watersheds to a target area identified by row-col indexes with labeling of nested watersheds
  IntegerMatrix watershed(flowdirs.nrow(), flowdirs.ncol());
  
  std::queue<cellz> expansion;
  
  for(int r = 0; r < target_rc.nrow(); r++){
    expansion.push(cellz(target_rc(r, 0)-1, target_rc(r, 1)-1, target_rc(r, 2)));
  }
  
  while(expansion.size()>0){
    cellz c = expansion.front();
    expansion.pop();
    
    for(int n = 1; n <= 8; n++){
      
      int nc = c.c+dx[n];
      int nr = c.r+dy[n];
      double label = 1;
      
      if(nested){
        label = c.z;
      }
      
      if(!(0<=nc && nc<flowdirs.ncol() && 0<=nr && nr<flowdirs.nrow()))
        continue;
      
      else if(flowdirs(nr, nc) == flowdir_nodata)
        continue;
      
      else if(watershed(nr, nc) == watershed_nodata && n == d8_inv[flowdirs(nr, nc)]){
        expansion.push(cellz(nr, nc, label));
        watershed(nr, nc) = label;
      }
    }
  }
  
  return watershed;
}



// // rc to i and vice versa works for column-major R matrixes now
// // [[Rcpp::export]]
// List rc_to_i(int r, int c, NumericMatrix m){
//   int ncol = m.ncol();
//   int nrow = m.nrow();
//   
//   int i = c*nrow+r;
//   
//   List result = List::create(_["i"] = i);
//   
//   return result;
//   
// }
// // [[Rcpp::export]]
// List i_to_rc(int i, NumericMatrix m){
//   int ncol = m.ncol();
//   int nrow = m.nrow();
//   
//   int r = i % nrow;
//   int c = i / nrow;
//   
//   List result = List::create(_["r"] = r, _["c"] = c);
//   
//   return result;
//   
// }


// row-col subscript to single index and vice versa
int rc_to_i(int row, int col, NumericMatrix m){
  
  assert(0<=col && col<m.ncol() && 0<=row && row<m.nrow());
  
  int i = col*m.nrow()+row;

  return i;
}

int i_to_r(int i, NumericMatrix m){
  int row = i % m.nrow();
  return row;
}

int i_to_c(int i, NumericMatrix m){
  int col = i / m.nrow();
  return col;
}

// [[Rcpp::export]]
NumericMatrix comp_breach_lindsay2016(NumericMatrix dem){
  //Complete breaching algorithm:
  //"Lindsay, J.B., 2016. Efficient hybrid breaching-filling sink removal methods for flow path enforcement in digital elevation models: Efficient Hybrid Sink Removal Methods for Flow Path Enforcement. Hydrological Processes 30, 846--857. doi:10.1002/hyp.10648"
  //As implemented in RichDEM

  int NO_BACK_LINK = numeric_limits<int>::max();

  int UNVISITED = 0;
  int VISITED = 1;
  int EDGE = 2;

  IntegerMatrix backlinks(dem.nrow(), dem.ncol());
  fill(backlinks.begin(), backlinks.end(), NO_BACK_LINK);
  IntegerMatrix visited(dem.nrow(), dem.ncol());
  LogicalMatrix pits(dem.nrow(), dem.ncol());

  int total_pits = 0;
  priority_queue<cellz, vector<cellz>, greater<cellz>> pq; //slightly different queue used in RichDEM using insertion order also

  //Seed the priority queue
  for(int r = 0; r < dem.nrow(); r++){
    for(int c = 0; c < dem.ncol(); c++){

      if(dem(r, c) == dem_nodata)
        continue;

      if(r==0 || c==0 || c == dem.ncol()-1 || r == dem.nrow()-1){
        pq.push(cellz(r, c, dem(r, c)));
        visited(r, c) = EDGE;
        continue;
      }

      double lowest_neighbour = numeric_limits<double>::max();
      for(int n = 1; n <= 8; n++){
        const int nc = c+dx[n];
        const int nr = r+dy[n];

        if(dem(nr, nc) == dem_nodata){
          pq.push(cellz(r, c, dem(r, c)));
          visited(r, c) = EDGE;
          goto nextcell;
        }

        lowest_neighbour = min(dem(nr, nc), lowest_neighbour);
      }

      if(dem(r, c) <= lowest_neighbour){
        dem(r, c) = lowest_neighbour;
        pits(r, c) = true;
        total_pits++;
      }

      nextcell:;
    }
  }

  while(!pq.empty()){

    const cellz c = pq.top();
    pq.pop();

    if(pits(c.r,c.c)){
      //Locate a cell that is lower than the pit cell, or an edge cell
      int cc = rc_to_i(c.r, c.c, dem);               //Current cell on the path
      int cc_r = i_to_r(cc, dem);
      int cc_c = i_to_c(cc, dem);
      double target_height = dem(c.r,c.c);                     //Depth to which the cell currently being considered should be carved

      //Trace path back to a cell low enough for the path to drain into it, or
      //to an edge of the DEM
      while(cc != NO_BACK_LINK && dem(cc_r, cc_c) >= target_height){
        dem(cc_r, cc_c) = target_height;
        cc = backlinks(cc_r, cc_c);
        cc_r = i_to_r(cc, dem);
        cc_c = i_to_c(cc, dem);
      }

      --total_pits;
      
      if(total_pits==0)
        break;
    }

    //Looks for neighbours which are either unvisited or pits
    for(int n = 1; n <= 8; n++){
      const int nc = c.c+dx[n];
      const int nr = c.r+dy[n];

      if(!(0<=nc && nc<dem.ncol() && 0<=nr && nr<dem.nrow()))
        continue;

      if(dem(nr, nc) == dem_nodata)
        continue;

      if(visited(nr, nc) != UNVISITED)
        continue;

      double my_e = dem(nr, nc);

      pq.push(cellz(nr, nc, my_e));
      visited(nr, nc) = VISITED;
      backlinks(nr, nc) = rc_to_i(c.r, c.c, dem);
    }
  }

  return dem;
}

