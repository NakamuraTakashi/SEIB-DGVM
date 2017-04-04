//Virtual Forest visualizer
//by Hisashi SATO (2015 Nov 13)

// +KI1.0 +KF50.0 +KFF50 <-- Example of the input for the option box for generating 50yrs frames
                       
        //Include files
#include "colors.inc"  
#include "shapes.inc"  
#include "textures.inc"
#include "Woods.inc"   
#include "stones.inc"  
#include "glass.inc"   
#include "metals.inc"  

        //Open data file
#fopen data_file "../code/output/forest.txt" read

//[For separate files]
//#fopen data_file concat(str(int((Flame-1)*clock+1),1,0),".txt") read
                       
        //Camera setting
// [1]SKEWEDVIEW
//camera      { location <100,120,30> look_at <55,60,0> sky<0,0,1> angle 60 } 

// [2]CLOSE VIEW
//camera      { location <15,40,15> look_at <15,15,0> sky<0,0,1> angle 60 }

// [3]OVERVIEW
// For Fukido
//camera      { location <362,200,400> look_at <362,-250,0> sky<0,0,1> angle 60 } 
//camera      { location <362,100,700> look_at <362,-322,0> sky<0,0,1> angle 60 } 
//camera      { location <362,-322,800> look_at <362,-322,0> sky<0,0,1> angle 60 } 
// Zoom
camera      { location <362,-200,20> look_at <362,-250,0> sky<0,0,1> angle 60 } 

// [4]UPVIEW
//camera      { location <0,0,1> look_at <50,60,15> sky<0,0,1> angle 60 }

// [5]ANIMATION (rotation view)
//camera      { location <15-cos(2*pi*clock)*30,15+sin(2*pi*clock)*30,20> look_at <15,15,0> sky<0,0,1> angle 60 }


//camera      { location <15,15,50> look_at <15,15,0> sky<0,0,1> angle 60 }

        //Scene setting
//light_source{ <0,100,500> color rgb 2.3 shadowless}
//light_source{ <724,100,500> color rgb 2.3 shadowless}

// Parallel  Lights
light_source  {
                         <1000,1000,4000>
                         color  rgb 2  
                         parallel
}
 

global_settings { ambient_light color rgb 1}
background  {color Silver }
object      {polygon{4,<0,0>,<724,0>,<724,-644>,<0,-644> 
                    texture{ T_Grnt15 finish{reflection 0.0} }
                    }
            }

//fog{
//    color White
//    fog_type 2
//    fog_alt 0.1
//    fog_offset 0.1
//    distance 0.5
//    rotate x*90
//    turbulence z*0.2
//    turb_depth 0.2
//}

//sky_sphere{
//    pigment{gradient z 
//    color_map{
//    [0.0 White*0.9 ]
//    [ 1.0 color rgb<0.3,0.4,1.2>]
//    }
//    }
//}

  
        //Specify data location
 #declare year_omit = clock -1;
 #while ( year_omit > 0 )

   #read  (data_file,tree_number)
   #while ( tree_number > 0 )
     #read(data_file, bole_x, bole_y, crown_x, crown_y, bole_h, foliage_h, bole_d, foliage_d, pft)
     #declare tree_number = tree_number - 1;
   #end

 #declare year_omit = year_omit - 1;
 #end
  
      
        //Text
  #declare Text = concat (str(clock,3,0)," year")
  
  text {
    ttf "timrom.ttf" Text 0.02, 0
    //ttf "timrom.ttf" Text 0.03, 0
    pigment { Blue }
    //rotate x*120
    rotate x*130
    //translate <26,45,35>    
    //translate <27,35.3,6>    
    translate <53,128,66.1>    
  }

        //Tree drawing
 #read  (data_file, tree_number)
 #while ( tree_number > 0 )
    #read(data_file, bole_x, bole_y, crown_x, crown_y, bole_h, foliage_h, bole_d, foliage_d, pft)
    
    //Bole
    #if (bole_h>0.0)
        object {cylinder { <bole_x, -bole_y, 0.0>, <bole_x, -bole_y, bole_h>, bole_d }
	        texture {T_Brass_5E}
	        finish  {reflection 0.0}
        }
    #end
    
    //Foliage
    #if (foliage_h>0.0)
        object {cylinder { <crown_x, -crown_y, bole_h>, <crown_x, -crown_y, bole_h+foliage_h>, foliage_d }
	     #switch(pft)
	        #case(1) //Tropical broad-leaved evergreen (1)
  	                texture{ pigment{ OrangeRed } }
	                finish{diffuse 1.0 crand 0.0 phong 1.0 reflection 0.1} 
	                #break
  	        #case(2) //Tropical broad-leaved evergreen (2)
  	                texture{ pigment{ ForestGreen } }
                    finish{diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0} 
  	                #break
  	        #case(3) //Tropical broad-leaved evergreen (3)
                    texture {T_Grnt23}
  	                finish{diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0} 
  	                #break 
  	        #case(4) //Tropical broad-leaved evergreen (4)
  	                texture {pigment { Red_Marble} }
                    finish{diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0} 
  	                #break
	        #case(4) //Tropical broad-leaved raingreen
	                texture {Yellow_Pine}
	                finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                #break
	        #case(3) //Temperate needle-leaved evergreen
	                texture {T_Yellow_Glass}   //texture {pigment { Jade} }
	                finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                #break
	        #case(7) //Temperate broad-leaved evergreen
	                texture {T_Stone18}
	                finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                #break
	        #case(2) //Temperate broad-leaved summergreen
	                texture {T_Winebottle_Glass}
	                finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                #break
	        #case(9) //Boreal needle-leaved evergreen
	                texture {pigment { Blue_Agate}}
	                finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                #break
	        #case(10) //Boreal needle-leaved summergreen
	                //texture {T_Vicksbottle_Glass}
              	                texture {T_Dark_Green_Glass}

	                texture {T_Stone41}
	                finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                #break
	        #else    //Boreal broad-leaved summergreen
	                texture {T_Dark_Green_Glass}
	                finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	     #end
        }
    #end
        
 #declare tree_number = tree_number - 1;
 #end
        
        //Close file
#fclose data_file
