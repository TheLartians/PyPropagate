
namespace lars{

	template <unsigned axis,class A,class R> void ring_derivative(const A &array,R & result,typename A::Scalar s){
          
          if(array.shape() != result.shape()) throw std::runtime_error("array shape does not match result");
           
	  auto ring = [=](typename A::Scalar v)->typename A::Scalar{
	    return fmod(fmod(v,s)+s,s) ;
	  };

          auto c = s/2;

	  result.for_all_indices([&](typename A::Index idx){ 
	    if(idx.template get<axis>() == 0){
	      auto ip1 = idx;
	      ip1.template set<axis>(idx.template get<axis>() + 1);
	      auto v0 = array(idx);
	      auto vp1 = ring(array(ip1) - v0 + c) - c;
	      result(idx) = vp1;
	    }
	    else if(idx.template get<axis>() == array.shape().template get<axis>() - 1){
	      auto im1 = idx;
	      im1.template set<axis>(idx.template get<axis>() - 1);
	      auto v0 = array(idx);
	      auto vm1 = ring(v0 - array(im1) + c) - c;
	      result(idx) = vm1;
	    }    
	    else{
	      auto im1 = idx, ip1 = idx;
	      im1.template set<axis>(idx.template get<axis>() - 1); 
	      ip1.template set<axis>(idx.template get<axis>() + 1);
	      auto v0 = array(idx);
	      auto vp1 = ring(array(ip1) - v0 + c);
	      auto vm1 = ring(array(im1) - v0 + c);
	      result(idx) = 0.5 * (vp1 - vm1);
	    } 
	  });
	  
	}

	template <unsigned axis,class A> typename A::Copy ring_derivative(const A &array,typename A::Scalar s){
          typename A::Copy result;
          result.resize(array.shape());
          ring_derivative<axis>(array,result,s);
          return result; 
        }
}

