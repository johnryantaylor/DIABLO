C******************************************************************************|
C these are the hooks required by the 'basic' version of diablo, even
C when not running the 'ensemble' flavor    
C
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine init_march

	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine init_adjoint

	return
	end 
		
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine save_stats_adjoint(final)
	logical final

	return
	end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine backwards_constants(TEMP1,TEMP2,TEMP3,TEMP4,TEMP5)
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5

	return
	end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine quasi_reversible_xz
	
	return
	end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine quasi_reversible_y

	return
	end


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine septa_implicit_solve_u
	
	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine quasi_rev_per(TEMP5,I,J,K)
c If direction = -1 then add to RHS of CN coefficients.
c If direction = +1 then add to implicit solve of 
	integer J,K,I
	real*8 TEMP5
	
	return
	end
