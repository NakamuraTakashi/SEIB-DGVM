//Virtual Forest visualizer
//by Hisashi SATO (2003 Sep 04)
     

        //Include files
#include "colors.inc"
#include "shapes.inc"
#include "textures.inc"
#include "Woods.inc"
#include "stones.inc"
#include "glass.inc"
#include "metals.inc"


        //Specizfy year of drawing
#declare year =   30;
                

        //Open data file
#fopen data_file "output_forest.txt" read


        //Scene setting
// [1]SKEWEDVIEW
//camera      { location <100,120,30> look_at <55,60,0> sky<0,0,1> angle 60 } 

// [2]CLOSE VIEW
//camera      { location <15,60,35> look_at <15,-10,0> sky<0,0,1> angle 60 }
//camera      { location <60,15,35> look_at <-10,15,0> sky<0,0,1> angle 60 }

// [3]OVERVIEW
//camera      { location <15,62,25> look_at <15,-25,0> sky<0,0,1> angle 60 } 
//camera      { location <15,75,60> look_at <15,-30,0> sky<0,0,1> angle 60 } 
//camera      { location <15,60,65> look_at <15,-10,0> sky<0,0,1> angle 60 } 
//camera      { location <50,15,30> look_at <0,15,0> sky<0,0,1> angle 60 } 
//camera      { location <25,75,30> look_at <25,20,0> sky<0,0,1> angle 60 } 
camera      { location <50,140,80> look_at <50,50,0> sky<0,0,1> angle 60 } 
//camera      { location <50,50,150> look_at <50,50,0> sky<0,0,1> angle 60 }
                                  
// [4]UPVIEW
//camera      { location <0,0,1> look_at <50,60,15> sky<0,0,1> angle 60 }
//camera      { location <10,10,0> look_at <10,10,50> sky<0,0,1> angle 60 }
//camera      { location <20,20,80> look_at <20,20,0> sky<0,0,1> angle 60 }


// [5]Side VIEW
//camera      { location <70,20,25> look_at <15,20,15> sky<0,0,1> angle 60 }
//camera      { location <20,120,50> look_at <20,15,30> sky<0,0,1> angle 60 }


light_source{ <50,160,100> color 2 shadowless}
light_source{ <50,-60,100> color 2 shadowless}

 

background  {color Silver }
object      {polygon{4,<0,0>,<100,0>,<100,100>,<0,100> 
                    texture{ T_Silver_1D finish{reflection 0.05} }
//                    texture {T_Grnt10}
                    }
            }




          //Specify data location
 #while ( year > 1 )

   #read  (data_file,tree_number)
   #while ( tree_number > 0 )
     #read(data_file, bole_x, bole_y, crown_x, crown_y, bole_h, foliage_h, bole_d, foliage_d, pft)
     #declare tree_number = tree_number - 1;
   #end

 #declare year = year - 1;
 #end
 
  
  
        //Tree drawing
#read(data_file,tree_number)
#while ( tree_number > 0 )
    #read(data_file, bole_x, bole_y, crown_x, crown_y, bole_h, foliage_h, bole_d, foliage_d, pft)
    
    //Bole
    #if (bole_h>0.0)
        object {cylinder { <bole_x, bole_y, 0.0>, <bole_x, bole_y, bole_h>, bole_d }
	        texture {T_Brass_5E}
	        finish  {reflection 0.0}
        }
    #end
    
    //Foliage
    #if (foliage_h>0.0)
        object {cylinder { <crown_x, crown_y, bole_h>, <crown_x, crown_y, bole_h+foliage_h>, foliage_d }
	     #switch(pft)
	        #case(1) //Tropical broad-leaved evergreen (1)
                        texture {T_Copper_4A}
	                finish{diffuse 1.0 crand 0.0 phong 1.0 reflection 0.1} 
	                #break
	        #case(2) //Tropical broad-leaved evergreen (2)
                        #texture {T_Grnt23}
                        texture {T_Dark_Green_Glass}
                        finish{diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0} 
	                #break
	        #case(3) //Tropical broad-leaved evergreen (3)
	                texture {Yellow_Pine}
	                finish{diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0} 
	                #break 
	        #case(4) //Tropical broad-leaved evergreen (4)
	                #texture {pigment { Red_Marble} }
                        texture {T_Winebottle_Glass}
                        finish{diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0} 
	                #break
	        #case(5)  //Tropical broad-leaved evergreen (Africa)
	                texture {T_Stone12}
                        finish{diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0} 
	                #break
	        #case(6)  //Tropical broad-leaved raingreen  (Africa)
	                texture {T_Ruby_Glass}
                        finish{diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0} 
	                #break
	        
	        #case(7) //Temperate needle-leaved evergreen
	                texture {T_Yellow_Glass}   //texture {pigment { Jade} }
	                finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                #break
	        #case(8) //Temperate broad-leaved evergreen
	                texture {T_Stone18}
	                finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                #break
	        #case(9)   //Temperate broad-leaved summergreen
	                texture {T_Winebottle_Glass}
	                finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                #break
	               	                
	        #case(10) //Boreal needle-leaved evergreen
	                texture {T_Winebottle_Glass}
	                finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                #break
	        #case(11) //Boreal needle-leaved summergreen
	                texture {T_Dark_Green_Glass}
	                finish{diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0} 
	                //finish {diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                #break
	        #else    //Boreal broad-leaved summergreen
	                texture {T_Vicksbottle_Glass}
	                finish{diffuse 0.5 crand 0.0 phong 1.0 reflection 0.0}
	                
	     #end
        }
    #end
        
    #declare tree_number = tree_number - 1;
#end

        
        
        //Close file
#fclose data_file
