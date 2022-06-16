    program SWE_1D_wetDry_nonhydro
      implicit none
      ! declare the variables -------------------

      ! integer
      integer i,j,m,step,step_n,total_steps,total_nodes,n_near,n_correct,i_correct,n_save,n_points
      integer, allocatable :: n_loc(:),loc_index(:,:),kin(:)

      ! real
      real*8, allocatable :: x(:),h(:,:),u(:,:),w(:,:),zeta(:,:),zeta_x(:,:),u_x(:,:),uh_x(:,:),x_c(:),q(:),q_x(:)
      real*8, allocatable :: weighting(:,:),lambda(:,:,:),zb(:),zb_x(:),u_bc(:),n_manning(:),x_point(:)
      real*8 g,t,ub,hb,rho,zetab, dx1, dx2

      ! string
      character*51 filename
      character*42 fileString

      ! parameter
      parameter(g=9.81d0)


      ! declare variables for initial condition

      real*8 dt, dx, h0, slope, h_min

      ! declare variables for wave maker

      real*8 amplitude, x0, c, wavemaker_eta, term1, term2
      real*8, allocatable :: wavemaker_hArr(:), wavemaker_uArr(:), ra(:)

      ! parameter control -------------------------

      m=3
      n_near=28
	  total_nodes = 1401
	  total_steps = 8000
	  n_correct = 5
	  dt = 0.001
	  dx = 0.002
	  h0 = 0.1
	  h_min = 0.001
	  rho = 1000
	  slope = 2.75

	  ! wave maker condition
	  amplitude = h0 * 0.12
	  x0 = -1.5
	  c = (g * (h0 + amplitude))**0.5

	  ! directory of the output data
      fileString = 'your directory'

	  !--------------------------------------------

      allocate(kin(1:total_nodes))
      allocate(x(1:total_nodes))
      allocate(x_c(1:total_nodes))
      allocate(h(1:total_nodes,0:2))
      allocate(u(1:total_nodes,0:2))
      allocate(w(1:total_nodes,0:2))
      allocate(zeta(1:total_nodes,0:2))
	  allocate(zeta_x(1:total_nodes,0:2))
      allocate(u_x(1:total_nodes,0:2))
      allocate(uh_x(1:total_nodes,0:2))
      allocate(zb(1:total_nodes))
	  allocate(zb_x(1:total_nodes))
	  allocate(q(1:total_nodes))
	  allocate(q_x(1:total_nodes))
      allocate(weighting(1:total_nodes,1:n_near*5))
      allocate(lambda(1:total_nodes,1:m,1:n_near*5))
      allocate(n_loc(1:total_nodes))
      allocate(loc_index(1:total_nodes,1:n_near*5))
      allocate(u_bc(0:total_steps))
      allocate(n_manning(1:total_nodes))
      allocate(ra(1:total_nodes))
      allocate(wavemaker_hArr(1:total_steps))
      allocate(wavemaker_uArr(1:total_steps))

      step = 0

      ! wave maker on the boundary (where x = 0)

      do j = 1,total_steps

        term1 = ( (3 * amplitude) / (4 * (h0**3) ) )**0.5
        term2 =  x0 + (c * j * dt)
        wavemaker_eta = amplitude * ( 1 / cosh(term1 * term2) )**2
        wavemaker_hArr(j) = wavemaker_eta + h0
        wavemaker_uArr(j) = (c * wavemaker_eta) / (wavemaker_eta + h0)

      end do

      ! initial condition set up

	  do j=1,total_nodes

		  ! x
		  if(j .eq. 1)then
            x(j) = 0
          else
            x(j) = x(j - 1) + dx
          end if

		  ! kin

		  if(j .eq. 1)then
		    kin(j) = 2
		  else if(j .eq. total_nodes) then
		    kin(j) = 1
		  else
		    kin(j) = 0
		  end if


		  ! zb and h --------------

		  if(j .le. ceiling(22.87 * h0 / dx) + 1)then

            h(j,0) = h0
            zb(j) = -h0

          else if(j .gt. ceiling(22.87 * h0 / dx) + 1 .and. j .le. ceiling(25.62 * h0 / dx) + 1)then

            h(j,0) = h(j - 1,0) - (dx / slope)
            zb(j) = -h(j,0)

          else

            h(j,0) = 0
            zb(j) = zb(j - 1) + (dx / slope)

          end if

          zeta(j,step) = h(j,0) + zb(j)
          u(j,0) = 0
          n_manning(j) = 0

	  end do


      ! calculate zb_x by definition instead of MLS, but you can toggle the comment
      ! of MLS_zb function if you want to obtain zb_x by MLS
	  do j = 1,total_nodes

        if(j .eq. total_nodes) then
            zb_x(j) = 0
        else
            zb_x(j) = (zb(j + 1) - zb(j)) / dx
        end if

      end do


	  ! start calculation -----------

      ! step 0 ----------------------
      step = 0
      step_n = 0

      call pre_MLS(total_nodes,m,n_near,x,weighting,lambda,n_loc,loc_index,g,h,dt,kin) !�p��U�I�������x�}�A�YLambda�x�}
	  !call MLS_zb(total_nodes,m,n_near,weighting,lambda,n_loc,loc_index,zb,zb_x)     !�p��zb�����ɼ�

      call MLS(total_nodes,m,n_near,weighting,lambda,n_loc,loc_index,zeta(:,step),h(:,step),u(:,step) &
              ,zeta_x(:,step),u_x(:,step),uh_x(:,step),zb,h_min,x,ra,zb_x,step,step_n)     !�p��h, u, uh�����ɼ�


	  call calculate_w(total_nodes,u(:,step),zeta_x(:,step),uh_x(:,step),zb_x,w(:,step))
	  call MLS_w(total_nodes,m,n_near,weighting,lambda,n_loc,loc_index,w(:,step),step_n)


      ! step 1 ----------------------
      step = 1 !�p���1�ɶ��B�H��A�Ystep>=1

      do step_n=1,total_steps  !�s�[�@�ӭp��index, step_n

	    ub = wavemaker_uArr(step_n)
	    hb = wavemaker_hArr(step_n)

        write(*,'(a5,i5)') 'step=',step_n

        call predict(total_nodes,step,dt,g,rho,n_manning,zb,zeta,zeta_x,h,u,u_x,uh_x,h_min)

	    call NF_BC(total_nodes,step,u,kin)

	    call WM_BC(total_nodes,step,h,u,kin,ub,hb)


        call MLS(total_nodes,m,n_near,weighting,lambda,n_loc,loc_index,zeta(:,step),h(:,step),u(:,step) &
                ,zeta_x(:,step),u_x(:,step),uh_x(:,step),zb,h_min,x,ra,zb_x, step_n,step_n)     !�p��h, u, uh�����ɼ�

        call calculate_w(total_nodes,u(:,step),zeta_x(:,step),uh_x(:,step),zb_x,w(:,step))


		call MLS_w(total_nodes,m,n_near,weighting,lambda,n_loc,loc_index,w(:,step),step_n)

		do i_correct=1,n_correct

	      call calculate_q(total_nodes,step,dt,h,w,q)

	      call MLS_q(total_nodes,m,n_near,weighting,lambda,n_loc,loc_index,q,q_x,step_n)     !�p��q�����ɼ�
          call correct(total_nodes,step,dt,g,rho,n_manning,zb,zb_x,zeta,zeta_x,h,u,u_x &
                ,uh_x,q,q_x,n_correct,i_correct,h_min,step_n)

	      call NF_BC(total_nodes,step,u,kin)

	      call WM_BC(total_nodes,step,h,u,kin,ub,hb)

          call MLS(total_nodes,m,n_near,weighting,lambda,n_loc,loc_index,zeta(:,step),h(:,step),u(:,step) &
                  ,zeta_x(:,step),u_x(:,step),uh_x(:,step),zb,h_min,x,ra,zb_x,step_n,step_n)     !�p��h, u, uh�����ɼ�

          call calculate_w(total_nodes,u(:,step),zeta_x(:,step),uh_x(:,step),zb_x,w(:,step))

		  call MLS_w(total_nodes,m,n_near,weighting,lambda,n_loc,loc_index,w(:,step),step_n)


		end do

2       do j=1,total_nodes
          zeta(j,0)=zeta(j,1)
          h(j,0)=h(j,1)
          u(j,0)=u(j,1)
          w(j,0)=w(j,1)
          zeta_x(j,0)=zeta_x(j,1)
          u_x(j,0)=u_x(j,1)
          uh_x(j,0)=uh_x(j,1)
        end do

!          ! data output
!          if( mod(step_n, 10) .eq. 0)then
!              write(filename,'(a42,i5.5,a4)') fileString, step_n, '.txt'
!
!              open(3,file=filename)
!
!                      do j = 1,total_nodes
!                          write(3,'(11f16.5)')x(j) , zeta(j, step)
!                      end do
!
!              close(3)
!          end if

      end do

    end


    subroutine pre_MLS(N,m,n_near,x,wt,lambda,n_loc,loc_index,g,h,dt,kin)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),kin(1:N)
      real*8 x(1:N),wt(1:N,1:n_near*5),lambda(1:N,1:m,1:n_near*5),d_near(1:n_near),p(1:m)
      real*8 r,x_loc,ro,g,dt,epsilon,h(1:N),rn
      real*8 a(1:n_near*5,1:m),at(1:m,1:n_near*5),ata(1:m,1:m),ata_inv(1:m,1:m),aa(1:15),xx(1:15)
	  aa(1)=3.09582566312239d0
	  aa(2)=41.7393749158713d0
	  aa(3)=-44.4172597684657d0
	  aa(4)=-38.5119790368627d0
	  aa(5)=-54.5886020275925d0
	  aa(6)=140.855346281139d0
	  aa(7)=40.5202718914226d0
	  aa(8)=45.0283055360917d0
	  aa(9)=46.883600490321d0
	  aa(10)=-186.882390021869d0
	  aa(11)=-14.0915137261954d0
	  aa(12)=-13.2222810955978d0
	  aa(13)=-20.6947245364581d0
	  aa(14)=-13.1601012713356d0
	  aa(15)=84.829579436323d0


  !�}�l��Collocation
      do j=1,N
  !���n_near�̪��I���Z��
        call KNN(N,x,d_near,n_near,j)
        ro=d_near(n_near)*1.01
!		if(kin(j) .eq. 0) then
!		  rn=(d_near(2)+d_near(3))/2
!		else
!		  rn=(d_near(2)+d_near(3))/3
!		end if
!		xx(1)=1
!		xx(2)=rn/h(j)/.3d0
!		xx(3)=(g*h(j))**.5/(rn/dt)/0.132060592153753d0
!		xx(4)=xx(2)**2
!		xx(5)=xx(2)*xx(3)
!		xx(6)=xx(3)**2
!		xx(7)=xx(2)**3
!		xx(8)=xx(2)**2*xx(3)
!		xx(9)=xx(2)*xx(3)**2
!		xx(10)=xx(3)**3
!		xx(11)=xx(2)**4
!		xx(12)=xx(2)**3*xx(3)
!		xx(13)=xx(2)**2*xx(3)**2
!		xx(14)=xx(2)*xx(3)**3
!		xx(15)=xx(3)**4
!		epsilon=.0d0
!		do i=1,15
!		  epsilon=epsilon+aa(i)*xx(i)
!		end do
!		epsilon=max(epsilon,0.d0)

        epsilon = 5

  !�j�I�íp���v��
        n_loc(j)=0
        do l=1,N !�v�@�j�I
          x_loc=x(l)-x(j) !�p��۹理��
          r=abs(x_loc)  !�p���focused node�������Z��
          if(r .ge. ro) goto 1  !�Y�Zfocused node�ӻ��A�]�n�ư�
          n_loc(j)=n_loc(j)+1
          loc_index(j,n_loc(j))=l !�o�̪�local_index�N�O�פ�̪�k�C�Y��l��node�b��j�ӧ����d�򤺡A��s����k
          wt(j,n_loc(j))=((1-r/ro)**epsilon)**.5 !�p���v���Y��
1       end do
      end do
!�غc�����x�}
      do j=1,N
  !�p��A�x�}�����e
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x(j)  !�p�⧽������
          call local_poly(x_loc,p,m) !�p��Local Polynomial
          do i=1,m
            a(k,i)=wt(j,k)*p(i)
          end do
        end do
  !����A�x�}����m�x�}
        do k=1,n_loc(j)
          do i=1,m
            at(i,k)=a(k,i)
          end do
        end do
  !����Psuudo�x�}
        do i=1,m
          do l=1,m
            ata(i,l)=0
            do k=1,n_loc(j)
              ata(i,l)=ata(i,l)+at(i,k)*a(k,l)    !atxa=at*a
            end do
          end do
        end do
  !����Psuudo�x�}���ϯx�}
        call matrix_inv(m,ata,ata_inv)
  !����lambda�x�}
        do i=1,m
          do k=1,n_loc(j)
            lambda(j,i,k)=0
            do l=1,m
              lambda(j,i,k)=lambda(j,i,k)+ata_inv(i,l)*at(l,k)
            end do
          end do
        end do
      end do
    end


! �o�ӰƵ{����X(XX)����K�̪��I���Z��
    subroutine KNN(N,x,d_near,n_near,star) ! K Nearest Neighbors
      implicit none
      integer N,i,j,k,n_near,star
      real*8 x(1:N),r,d_near(1:n_near)
  !  �@�}�l�̪�Z�����]1.d99
      do i=1,n_near
        d_near(i)=1.d99
      end do
  !  �}�l�j�M
      do j=1,N
        r=abs(x(j)-x(star))
        do i=1,n_near
          if(r .lt. d_near(i)) then
            do k=n_near-1,i,-1
              d_near(k+1)=d_near(k)
            end do
            d_near(i)=r
            goto 1
          end if
        end do
1     end do
    end




    subroutine local_poly(x,p,m)
      implicit none
      integer m
      real*8 x,p(1:m)
      p(1)=1.0d0
      p(2)=x
      p(3)=x**2/2
      !p(4)=x**3/6
    end




    subroutine matrix_inv(m,a,a_inv)
      implicit none
      integer i,j,k,m
      real*8 a(1:m,1:m),a_inv(1:m,1:m),z
      do i=1,m
        do j=1,m
          if(i .eq. j) then
            a_inv(i,j)=1.d0
          else
            a_inv(i,j)=0
          end if
        end do
      end do
      do i=1,m-1
        z=a(i,i)
        do j=i,m
          a(i,j)=a(i,j)/z
        end do
        do j=1,m
          a_inv(i,j)=a_inv(i,j)/z
        end do
        do k=i+1,m
          z=a(k,i)
          do j=1,m
            a(k,j)=a(k,j)-z*a(i,j)
            a_inv(k,j)=a_inv(k,j)-z*a_inv(i,j)
          end do
        end do
      end do
      do j=1,m
        a_inv(m,j)=a_inv(m,j)/a(m,m)
      end do
      do i=m-1,1,-1
        do j=m,i+1,-1
          do k=1,m
            a_inv(i,k)=a_inv(i,k)-a_inv(j,k)*a(i,j)
          end do
        end do
        do k=1,m
          a_inv(i,k)=a_inv(i,k)/a(i,i)
        end do
      end do
    end



    subroutine MLS_zb(N,m,n_near,wt,lambda,n_loc,loc_index,zb,zb_x)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5)
      real*8 x(1:N),zb(1:N),zb_x(1:N),new_zb(1:N),lambda(1:N,1:m,1:n_near*5),wt(1:N,1:n_near*5)
      real*8 alpha(1:m),beta(1:n_near*5)
  !�}�l�p��
    ! zb������
      do j=1,N
      !�p��beta_zb
        do k=1,n_loc(j)
          beta(k)=wt(j,k)*zb(loc_index(j,k))
        end do
      !�p��alpha_zb
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !��X zb,zb_x
        !new_zb(j)=alpha(1)
        zb_x(j)=alpha(2)
      end do
!  !�N��
!      do j=1,N
!        zb(j)=new_zb(j)
!      end do
    end



    subroutine MLS(N,m,n_near,wt,lambda,n_loc,loc_index,zeta,h,u,zeta_x,u_x,uh_x,zb,h_min,x,ra,zb_x,step,step_n)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),step,step_n
      real*8 zeta(1:N),h(1:N),u(1:N),w(1:N),new_zeta(1:N),new_u(1:N) &
            ,zeta_x(1:N),u_x(1:N),uh_x(1:N),lambda(1:N,1:m,1:n_near*5),wt(1:N,1:n_near*5)
      real*8 alpha(1:m),beta(1:n_near*5),zb(1:N)
      ! parameter for wet dry
      real*8 wet(1:N), shore(1:N), h_min, r, ra(1:N), r_nearest_dry(1:N)
      real*8 x(1:N), x_loc, wet_neighbors(1:N), ww
      real*8 theta(1:N), w_e(1:n_near*5), lambda_e(1:m,1:n_near*5)
      real*8 p(1:m) ,a(1:n_near*5,1:m), at(1:m,1:n_near*5), ata(1:m,1:m), ata_inv(1:m,1:m)
      real*8 zb_x(1:N)

      ! output wetDry condition, please change the directory to the place you want to output
      ! Notice that the length of the characters should also be changed
      character*44 fileString
      character*53 filename

      fileString = 'your directory'

   !�ʲ��P�_����
      do j=1,N

        if(h(j) .le. h_min) then
          wet(j)=0
        else
          wet(j)=1
        end if

      end do

   !��X��F���I�����I�A�íp��Z�̪��I���h��
      do j=1,N

        shore(j)=0  !�����]��0
        r_nearest_dry(j)=1.d99 !�����]���ܤj

        do k=1,n_loc(j)

          x_loc=x(loc_index(j,k))-x(j)
          r = x_loc !(x_loc**2)**.5d0

          if(r .lt. ra(j)*4.d0/3.d0) then                                      !��F���I�����A
            if(wet(j) .eq. 0 .and. wet(loc_index(j,k)) .eq. 1) then            !�����I�A�B�o�@�I�����O���I
              if( zb(j)+h(j) .lt. zb(loc_index(j,k))+h(loc_index(j,k)) ) then  !�ӥ������I�����{���F���I�������٧C
                shore(j)=1                                                      !�h���O��F���I�����I�A�n�����I�B�z
              end if
            end if
          end if

          if(r .lt. r_nearest_dry(j) .and. wet(loc_index(j,k)) .eq. 0) then    !��F���I�����A�O���I�B�Z�����̪��I�٤p
            r_nearest_dry(j)=r                                                 !�h��s��̪��I���Z��
          end if

        end do

      end do

   !���s�w�q���I
      do j=1,N

        wet(j)=wet(j)+shore(j) !���O���I�A�h����1�F�Y�O���I�A���n�����I�B�z�A�h0+1=1�F���I�Z���I�ܻ��A�h��0

      end do


   !�p��ƶq�Τ�v
      do j=1,N
        wet_neighbors(j)=0
        do k=1,n_loc(j)
          if(wet(loc_index(j,k)) .eq. 1) then
            wet_neighbors(j)=wet_neighbors(j)+1
          end if
        end do
        !theta(j)=wet_neighbors(j)/n_loc(j) ! theta �N��P�����I����v�A�Y��1�h��ܩP�򳣬O���I�A�Y��0�h��ܩP�򳣬O���I
        theta(j)=1.d0*wet_neighbors(j)/n_loc(j)
      end do



   !�����I�w�qzeta�ȡC���I��zeta�ä��s�b�A�γ̪����I��zeta�ȷ���zeta�C�Y�̪����I����ӥH�W�A�h�������C�Y�d�򤺧����L���I�A�h�Ω��ɰ��{��zeta�ȡC
      do j=1,N

        if(theta(j) .eq. 0) then
          zeta(j) = 1.99  !�Y�d�򤺧����L���I�A�h�]zeta=1.99�A�b�p��zeta���ɼƮɡA�o�����I�n�ư��A�Y�n���s�غc�����x�}
          goto 11
        end if

        if(wet(j) .eq. 1) then
          zeta(j) = zb(j)+h(j)  !�Y�����O���I�A�h��zeta=zb+h
          !zeta(j) = zeta(j)
          goto 11
        end if

        zeta(j) = 0  !�ѤU���N�O���I
        ww = 0
        do k=1,n_loc(j)

          if(wet(loc_index(j,k)) .eq. 1) then
             ww = ww + wt(j,k)**2
             zeta(j) =  zeta(j) + ( h(loc_index(j,k))+zb(loc_index(j,k)) ) * wt(j,k)**2
          end if

        end do

         zeta(j) =  zeta(j)/ww

11    end do

!      ! wet dry condition data output and the length of the file name should also be changed
!      if (step_n .gt. 0) then
!
!          write(filename,'(a44,i5.5,a4)') fileString, step_n, '.txt'
!
!          open(3,file=filename)
!
!                  do j = 1,N
!                      write(3,'(11f16.5)')x(j) , wet(j)
!                  end do
!
!          close(3)
!
!
!      end if


  !�}�l�p��
    ! zeta������
      do j=1,N

      !�p��beta_zeta
        do k=1,n_loc(j)

          if(wet(j) .eq. 0 .or. theta(j) .lt. .15d0) then !�Y�o�@�I�O���I�A�νd�����I�ƶq�ܤ֡A�h�]�����y�ʡC
             new_zeta(j) = zb(j)
             zeta_x(j) = zb_x(j)
            goto 21
          end if

          if(zeta(loc_index(j,k)) .gt. 1.98) then
            goto 22 !�Y�d�򤺦�eta=1.d99���I�A�h�n���s�غc�����x�}
          end if

          beta(k)=wt(j,k)*zeta(loc_index(j,k))

        end do

      !�p��alpha_zeta
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !��X zeta, zeta_x
        new_zeta(j)=alpha(1)
        zeta_x(j)=alpha(2)
        goto 21

      !���s�غc�����x�}
      !�p��A�x�}�����e
22      do k=1,n_loc(j)
          if(zeta(loc_index(j,k)) .gt. 1.98) then
            w_e(k)=0 !�Yeta=1.d99�A�h�ư��A�]�v����0
          else
            w_e(k)=wt(j,k) !���M�N�ӭ쥻���v��
          end if
          x_loc=x(loc_index(j,k))-x(j)
          call local_poly(x_loc,p,m) !�p��Local Polynomial
          do i=1,m
            a(k,i)=w_e(k)*p(i)  ! ��eta�ɡA�u��eta������1.d99���I�Ӻ�C
          end do
        end do
      !����A�x�}����m�x�}
        do k=1,n_loc(j)
          do i=1,m
            at(i,k)=a(k,i)
          end do
        end do
      !����Psuudo�x�}
        do i=1,m
          do l=1,m
            ata(i,l)=0
            do k=1,n_loc(j)
              ata(i,l)=ata(i,l)+at(i,k)*a(k,l)    !atxa=at*a
            end do
          end do
        end do
        call matrix_inv(m,ata,ata_inv)
      !����lambda�x�}
        do i=1,m
          do l=1,n_loc(j)
            lambda_e(i,l)=0
            do k=1,m
              lambda_e(i,l)=lambda_e(i,l)+ata_inv(i,k)*at(k,l)
            end do
          end do
        end do
      !�p��beta_zeta
        do k=1,n_loc(j)
          beta(k)=w_e(k)*zeta(loc_index(j,k))
        end do
      !�p��alpha_zeta
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda_e(i,k)*beta(k)
          end do
        end do
      !��X zeta, zeta_x
        new_zeta(j)=alpha(1)
        zeta_x(j)=alpha(2)

21      end do

    ! u uh ������
      do j=1,N

      !�Y�o�@�I�O���I�A�νd�����I�ƶq�ܤ֡A�h�]�����y��
        if(wet(j) .eq. 0 .or. theta(j) .lt. .15d0) then
          new_u(j) = 0
          u_x(j) = 0
          uh_x(j) = 0
          goto 31
        end if

      !�p��beta_u
        do k=1,n_loc(j)

            beta(k)=wt(j,k)*u(loc_index(j,k))

        end do

      !�p��alpha_u
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !��X u,u_x
        new_u(j)=alpha(1)
        u_x(j)=alpha(2)

      !�p��beta_uh
        do k=1,n_loc(j)

            beta(k)=wt(j,k)*u(loc_index(j,k))*h(loc_index(j,k))

        end do
      !�p��alpha_uh
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !��X uh_x
        uh_x(j)=alpha(2)

31	  end do


   !�N��
      do j=1,N
        zeta(j)=new_zeta(j)
        u(j)=new_u(j)
		h(j)=zeta(j)-zb(j)

		if(h(j) .lt. 0) then
          h(j)=0
          zeta(j) = zb(j)
        end if

        if(h(j) .le. h_min) then
          u(j)=0
        end if

      end do

    end


    subroutine calculate_w(total_nodes,u,zeta_x,uh_x,zb_x,w)
      implicit none
      integer j,total_nodes
      real*8 u(1:total_nodes),zeta_x(1:total_nodes),uh_x(1:total_nodes),zb_x(1:total_nodes),w(1:total_nodes)

      do j=1,total_nodes
        w(j)=.5d0*(u(j)*(zeta_x(j)+zb_x(j))-uh_x(j))
      end do

	end



    subroutine predict(total_nodes,step,dt,g,rho,n,zb,zeta,zeta_x,h,u,u_x,uh_x,h_min)
      implicit none
      integer j,total_nodes,step
      real*8 zb(1:total_nodes),zeta(1:total_nodes,0:2),h(1:total_nodes,0:2),u(1:total_nodes,0:2) &
            ,zeta_x(1:total_nodes,0:2),u_x(1:total_nodes,0:2),uh_x(1:total_nodes,0:2) &
            ,dt,g,n(1:total_nodes),rho
      real*8 h_min
      do j=1,total_nodes

        zeta(j,step)=zeta(j,step-1)-dt*uh_x(j,step-1)

        h(j,step)=zeta(j,step)-zb(j)

        if(h(j,step) .lt. 0) then
            h(j,step)=0
        end if

        if(h(j,step) .le. h_min) then
            u(j,step)=0
            goto 1
        end if

        u(j,step)=u(j,step-1)-dt*(g*zeta_x(j,step-1) &
		                         +u(j,step-1)*u_x(j,step-1))
								 !+n(j)**2*g*u(j,step-1)*abs(u(j,step-1))/h(j,step-1)**(4.d0/3))

1     end do
    end




    subroutine correct(total_nodes,step,dt,g,rho,n,zb,zb_x,zeta,zeta_x,h,u,u_x,uh_x,q,q_x &
        ,n_correct,i_correct, h_min,step_n)
      implicit none
      integer j,total_nodes,step,n_correct,i_correct, step_n
      real*8 zb(1:total_nodes),zb_x(1:total_nodes),zeta(1:total_nodes,0:2),h(1:total_nodes,0:2) &
	        ,u(1:total_nodes,0:2),q(1:total_nodes),zeta_x(1:total_nodes,0:2),q_x(1:total_nodes) &
            ,h_x(1:total_nodes,0:2),u_x(1:total_nodes,0:2),uh_x(1:total_nodes,0:2) &
            ,dt,g,s1,s2,s3,s4,s5,rr,n(1:total_nodes),rho,u_mean,h_mean
      real*8 h_min

      rr=min(1.d0*i_correct/(n_correct-1),1.d0)
      rr=1-abs(cos(3.14159265358979*rr/2))**2.d0

      do j=1,total_nodes

        zeta(j,step)=zeta(j,step-1)-dt*(uh_x(j,step-1)+uh_x(j,step))/2
        h(j,step)=zeta(j,step)-zb(j)

        if(h(j,step) .lt. 0) then
            h(j,step) = 0
            zeta(j,step) = 0
        end if

        if(h(j,step) .le. h_min) then
            u(j,step)=0
            goto 1
        end if

	    h_mean=(h(j,step)+h(j,step-1))/2
		u_mean=(u(j,step)+u(j,step-1))/2
		s1=g*(zeta_x(j,step-1)+zeta_x(j,step))/2
		s2=(u(j,step-1)*u_x(j,step-1)+u(j,step)*u_x(j,step))/2
        s3=.5d0*q_x(j)*rr
		s4=.5d0*q(j)*rr*(zeta_x(j,step-1)+zeta_x(j,step)+2*zb_x(j))/(2*h_mean)
		s5=n(j)**2*g*u_mean*abs(u_mean)/h_mean**(4.d0/3)
		u(j,step)=u(j,step-1)-dt*(s1+s2+s3+s4)

1      end do
    end


    subroutine NF_BC(total_nodes,step,u,kin)
      implicit none
      integer j,total_nodes,step,kin(1:total_nodes)
      real*8 h(1:total_nodes,0:2),u(1:total_nodes,0:2)
      do j=1,total_nodes
        if(kin(j) .eq. 1) then
          u(j,step)=0
        end if
      end do
    end




    subroutine calculate_q(total_nodes,step,dt,h,w,q)
      implicit none
      integer j,total_nodes,step
      real*8 h(1:total_nodes,0:2),w(1:total_nodes,0:2),q(1:total_nodes),dt,w_t
      do j=1,total_nodes
	    w_t=(w(j,step)-w(j,step-1))/dt
        q(j)=(h(j,step)+h(j,step-1))/2*w_t
      end do
    end




    subroutine MLS_q(N,m,n_near,wt,lambda,n_loc,loc_index,q,q_x,step_n)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),n_smooth,step_n
      real*8 q(1:N),q_x(1:N),lambda(1:N,1:m,1:n_near*5),wt(1:N,1:n_near*5)
      real*8 alpha(1:m),beta(1:n_near*5),new_q(1:N),aa

      n_smooth = 8

	  do i=1,n_smooth
	    do j=1,N
	      aa=0
		  new_q(j)=0
          do k=1,n_loc(j)
		    aa=aa+wt(j,k)**2
		    new_q(j)=new_q(j)+wt(j,k)**2*q(loc_index(j,k))
		  end do
		  new_q(j)=new_q(j)/aa
	    end do
	    do j=1,N
	      q(j)=new_q(j)
	    end do
	  end do
  !�}�l�p��
    ! q������
      do j=1,N
      !�p��beta_q
        do k=1,n_loc(j)
          beta(k)=wt(j,k)*q(loc_index(j,k))
        end do
      !�p��alpha_q
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !��X q_x
        new_q(j)=alpha(1)
        q_x(j)=alpha(2)
      end do
  !�N��
      do j=1,N
        q(j)=new_q(j)
      end do
    end



    subroutine WM_BC(total_nodes,step,h,u,kin,ub,hb)
      implicit none
      integer j,total_nodes,step,kin(1:total_nodes)
      real*8 h(1:total_nodes,0:2),u(1:total_nodes,0:2),ub,hb
      do j=1,total_nodes
        if(kin(j) .eq. 2) then
          u(j,step)=ub
          h(j,step)=hb
          !u(j,step)=0
        end if
      end do
    end


    subroutine MLS_w(N,m,n_near,wt,lambda,n_loc,loc_index,w,step_n)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),n_smooth,step_n
      real*8 w(1:N),lambda(1:N,1:m,1:n_near*5),wt(1:N,1:n_near*5)
      real*8 alpha(1:m),beta(1:n_near*5),new_w(1:N),aa

      n_smooth = 8

	  do i=1,n_smooth
	    do j=1,N
	      aa=0
		  new_w(j)=0
          do k=1,n_loc(j)
		    aa=aa+wt(j,k)**2
		    new_w(j)=new_w(j)+wt(j,k)**2*w(loc_index(j,k))
		  end do
		  new_w(j)=new_w(j)/aa
	    end do
	    do j=1,N
	      w(j)=new_w(j)
	    end do
	  end do

  !�}�l�p��
      do j=1,N
      !�p��beta_w
        do k=1,n_loc(j)
          beta(k)=wt(j,k)*w(loc_index(j,k))
        end do
      !�p��alpha_w
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !��X smooth_w
        new_w(j)=alpha(1)
      end do
  !�N��
      do j=1,N
        w(j)=new_w(j)
      end do
    end


