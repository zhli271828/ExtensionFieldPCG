#pragma once
#include <cryptoTools/Common/block.h>
#include <emmintrin.h>  
#include <immintrin.h>  
#include <bitset>


namespace osuCrypto{
   template<typename T>
	void inv_foleageFft(
		span<T> coeffs,
		const size_t num_vars,
		const size_t num_coeffs)
	{
		

		if (num_vars > 1)
		{
			
			inv_foleageFft(
				coeffs,
				num_vars - 1,
				num_coeffs / 3);

			
			inv_foleageFft(
				coeffs.subspan(num_coeffs),
				num_vars - 1,
				num_coeffs / 3);

			inv_foleageFft(
				coeffs.subspan(2 * num_coeffs),
				num_vars - 1,
				num_coeffs / 3);
		}

		// temp variables to store intermediate values
		T tL, tM;
		T mult, xor_h, xor_l;

		T* coeffsL = coeffs.data() + 0;
		T* coeffsM = coeffs.data() + num_coeffs;
		T* coeffsR = coeffs.data() + 2 * num_coeffs;

		T mask_h; // 1010101001010
		T mask_l; // 10101010100101
		setBytes(mask_h, 0b10101010);
		setBytes(mask_l, 0b01010101);

		for (size_t j = 0; j < num_coeffs; j++)
		{
			xor_h = (coeffsM[j] ^ coeffsR[j]) & mask_h;
			xor_l = (coeffsM[j] ^ coeffsR[j]) & mask_l;

			
			mult = ((xor_h >> 1) ^ xor_l) | (xor_l << 1);

		
			tL = coeffsL[j] ^ coeffsM[j] ^ coeffsR[j];

			
			tM = coeffsL[j] ^ coeffsR[j] ^ mult;

			coeffsM[j] = coeffsL[j] ^ coeffsM[j] ^ mult;

			coeffsL[j] = tL;
			coeffsR[j] = tM;
		}
	}



block interleave(uint64_t a, uint64_t b)
{
    const uint64_t EVEN_MASK = 0xAAAAAAAAAAAAAAAAULL;
    const uint64_t ODD_MASK  = 0x5555555555555555ULL;
    uint64_t even_part_low = _pdep_u64(b, EVEN_MASK);
    uint64_t odd_part_low  = _pdep_u64(a, ODD_MASK);
    uint64_t low64 = even_part_low | odd_part_low;
    uint64_t  high64 = _pdep_u64(b >> 32, EVEN_MASK) | _pdep_u64(a >> 32, ODD_MASK);
   // block temp(high64, low64);
   return block(high64, low64);
}

uint64_t odd_bit_extraction(block x)
{
   const uint64_t odd_mask = 0x5555555555555555ULL;
   uint64_t x_low=_mm_extract_epi64(x,0);
   uint64_t x_high=_mm_extract_epi64(x,1);
   uint64_t tempA=_pext_u64(x_low, odd_mask);
   uint64_t tempB=_pext_u64(x_high, odd_mask);
   return (tempB<<32) | tempA;

}

uint64_t even_bit_extraction(block x)
{
  const uint64_t even_mask = 0xAAAAAAAAAAAAAAAAULL;
  uint64_t x_low=_mm_extract_epi64(x,0);
  uint64_t x_high=_mm_extract_epi64(x,1);
  uint64_t tempA=_pext_u64(x_low, even_mask);
  uint64_t tempB=_pext_u64(x_high, even_mask);
  return (tempB<<32) | tempA;

}

block F4_reduce(block odd, block even)
{
    block mod(0b1010, 0b1); //(\alpha h2(x), h1(x)), odd=(f1(x)x^64+f2(x)), even=\alpha(g1(x)x^{64}+g2(x))
    block t1= _mm_clmulepi64_si128(odd, mod, 0x01); //t1(x)=h1(x)f1(x)
    block t2= _mm_clmulepi64_si128(odd, mod, 0x11);// t2(x)=h2(x)f1(x)
    block t3=_mm_clmulepi64_si128(even, mod, 0x01); //t3(x)=g1(x)h1(x)
    block t4=_mm_clmulepi64_si128(even, mod, 0x11);//t4(x)=g1(x)h2(x)
    //std::cout<<"t1="<<std::bitset<64>(_mm_extract_epi64(t1,0))<<std::endl;
    block reduced_odd=t1^t4;   
    block reduced_even=t2^t3^t4;
    //do it again to finally reduced the degree to less than x^64.
    t1=_mm_clmulepi64_si128(reduced_odd, mod, 0x01); //t1(x)=h1(x)f1(x)
    t2=_mm_clmulepi64_si128(reduced_odd, mod, 0x11);// t2(x)=h2(x)f1(x)
    t3=_mm_clmulepi64_si128(reduced_even, mod, 0x01); //t3(x)=g1(x)h1(x)
    t4=_mm_clmulepi64_si128(reduced_even, mod, 0x11);//t4(x)=g1(x)h2(x)
    block final_reduced_odd=t1^t4;
    block final_reduced_even=t2^t3^t4;//final_reduced should only contain degree less than x^64.
    uint64_t odd_64=_mm_extract_epi64(odd^final_reduced_odd^reduced_odd,0);
    uint64_t even_64=_mm_extract_epi64(even^final_reduced_even^reduced_even,0); 
    return interleave(odd_64, even_64);

}


block F4_mul(block x, block y)
{
   block prod_odd, prod_even;
   uint64_t x_odd=odd_bit_extraction(x);
   uint64_t x_even=even_bit_extraction(x);
   uint64_t y_odd=odd_bit_extraction(y);
   uint64_t y_even=even_bit_extraction(y);
   block temp_x=_mm_set_epi64x(x_even, x_odd);
   block temp_y=_mm_set_epi64x(y_even, y_odd);
   block t1 = _mm_clmulepi64_si128(temp_x, temp_y, 0x00);
   block t2 = _mm_clmulepi64_si128(temp_x, temp_y, 0x10);
   block t3 = _mm_clmulepi64_si128(temp_x, temp_y, 0x01);
   block t4 = _mm_clmulepi64_si128(temp_x, temp_y, 0x11);
   prod_odd=(t1 ^ t4);
   prod_even=((t2 ^ t3)^ t4);
   return F4_reduce(prod_odd, prod_even);
}
}





