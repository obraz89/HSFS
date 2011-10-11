#include "Elems.h"
#include <vector>
#include <iostream>
#ifndef __STAB_ELEMS
#define __STAB_ELEMS
// instability wave characteristics 
// 
	struct t_WaveChars{
		t_CompVal a,b,w;
		t_CompVal vga, vgb;	// group velocity
		t_CompVal resid;
		inline static t_WaveChars find_max_instab(const std::vector<t_WaveChars>& vec){
			if(vec.size()>0){
				const t_WaveChars* pmax = &vec[0];
				for (int k=0; k<vec.size(); k++){
					if (vec[k].w.imag()>pmax->w.imag()){
						pmax = &vec[k];
					}
				}
				return *pmax;
			}else{
				return t_WaveChars();
			}
		};
		void print(){
			std::cout<<"a:"<<this->a<<std::endl<<"b:"<<this->b<<std::endl<<"w:"<<this->w<<std::endl;
		}
	};

#endif // __STAB_ELEMS