    program SWE1D_wetDry_nonhydro
      implicit none
      integer i,j,m,step,step_n,total_steps,total_nodes,n_near,n_flux_BC,n_eta_BC,points,n_step_save
      integer count
      real*8 x_point(1:6)
      real*8, allocatable :: x(:),h(:,:),uh(:,:),h_x(:,:),uh_x(:,:),x_c(:) &
                            ,u2h_x(:,:),n_x(:),n_y(:),zb(:),zb_x(:),eta(:,:),eta1(:) &
                            ,flux_x(:,:),flux_y(:,:),flux1_x(:),flux1_y(:),coriolis(:) &
                            ,uh_xx(:,:),ra(:),wavemaker_uhArr(:),wavemaker_hArr(:)
      integer, allocatable :: kin(:),cor(:),region(:),nearest_b(:),nearest_i(:)
      real*8, allocatable :: weighting(:,:),lambda(:,:,:)
      integer, allocatable :: n_loc(:),loc_index(:,:)
      real*8 dt,dx,g,Cz,h_min,t_ramp,t,u_max,volume,h0,slope
      real*8 amplitude,x0,c,wavemaker_eta,wavemaker_u,term1,term2,wavemaker_uh_bc,wavemaker_h_bc
      real*8, allocatable :: compareWithexperiment(:),animation(:)
      real*8, allocatable :: w_velocity(:,:),q(:),q_x(:)
      integer n_correct,i_correct
      integer gaugePosition
      character*22 filename
      parameter(g=9.81d0)


	  ! parameter control -------------------------

      m=3
      n_near=10
	  total_nodes = 3000     ! 1376 2902
	  total_steps = 1000
	  n_correct = 1
	  n_eta_BC = 0
	  n_flux_BC = 0
	  t_ramp = 0
	  dt = 0.0005
	  dx = 0.002
	  h0 = 0.218
	  h_min = 0.001
	  slope = 2.75
	  !Cz = 55.56
	  Cz = 0
	  u_max = 25
	  points = 0
	  n_step_save = 250
	  gaugePosition = 1635

	  ! wave maker condition
	  amplitude = h0 * 0.12
	  x0 = -1.5
	  c = (g * (h0 + amplitude))**0.5

	  !--------------------------------------------

      allocate(kin(1:total_nodes))
      allocate(cor(1:total_nodes))
      allocate(region(1:total_nodes))
      allocate(x(1:total_nodes))
      allocate(x_c(1:total_nodes))
      allocate(h(1:total_nodes,0:2))
      allocate(uh(1:total_nodes,0:2))
      allocate(h_x(1:total_nodes,0:2))
      allocate(uh_x(1:total_nodes,0:2))
      allocate(u2h_x(1:total_nodes,0:2))
      allocate(n_x(1:total_nodes))
      allocate(zb(1:total_nodes))
      allocate(zb_x(1:total_nodes))

      ! non-hydrostatic
      allocate(w_velocity(1:total_nodes,0:2))
      allocate(q(1:total_nodes))
      allocate(q_x(1:total_nodes))

      allocate(eta(1001:1000+n_eta_BC,0:total_steps))
      allocate(eta1(1001:1000+n_eta_BC))
      allocate(flux_x(5001:5000+n_flux_BC,0:total_steps))
      allocate(flux1_x(5001:5000+n_flux_BC))
      allocate(weighting(1:total_nodes,1:n_near*5))
      allocate(lambda(1:total_nodes,1:m,1:n_near*5))
      allocate(n_loc(1:total_nodes))
      allocate(loc_index(1:total_nodes,1:n_near*5))
      allocate(coriolis(1:total_nodes))
      allocate(uh_xx(1:total_nodes,0:2))
      allocate(ra(1:total_nodes))
      allocate(nearest_b(1:total_nodes))
      allocate(nearest_i(1:total_nodes))
      allocate(wavemaker_uhArr(0:total_steps))
      allocate(wavemaker_hArr(0:total_steps))
      allocate(compareWithexperiment( 1:(total_steps) ) )
      allocate(animation(1:(total_steps * total_nodes)))

      ! wave maker on the boundary (where x = 0)

      do j = 0,total_steps

        term1 = ( (3 * amplitude) / (4 * (h0**3) ) )**0.5
        term2 =  x0 + (c * j * dt)
        wavemaker_eta = amplitude * ( 1 / cosh(term1 * term2) )**2
        wavemaker_u = (c * wavemaker_eta) / (wavemaker_eta + h0)
        wavemaker_uhArr(j) = wavemaker_u * (wavemaker_eta + h0)
        wavemaker_hArr(j) = wavemaker_eta + h0

      end do

	  ! initial data

	  do j=1,total_nodes
	      cor(j) = 0
		  region(j) = 0
		  coriolis(j) = 0
		  n_x(j) = 0

		  ! x
		  if(j .eq. 1)then
            x(j) = 0
          else
            x(j) = x(j - 1) + dx
          end if

		  ! kin
		  if(j .eq. 1)then
		    kin(j) = 2
		  else if(j .eq. total_nodes-1) then
		    kin(j) = 1
		  else if(j .eq. total_nodes) then
		    kin(j) = -1
		  else
		    kin(j) = 0
		  end if

		  ! zb and h
		  if(j .le. ceiling(22.87 * h0 / dx) + 1)then  ! 8.87 22.87

            h(j,0) = h0
            zb(j) = -h(j,0)

          else if(j .gt. ceiling(22.87 * h0 / dx) + 1 .and. j .le. ceiling(25.62 * h0 / dx) + 1)then  ! 11.62 25.62

            h(j,0) = h(j - 1,0) - (dx / slope)
            !zb(j) = zb(j - 1) + (dx / slope)
            zb(j) = -h(j,0)

          else

            h(j,0) = 0
            zb(j) = zb(j - 1) + (dx / slope)

          end if

		  uh(j,0) = 0

	  end do

	  ! start calculation -----------
      call pre_MLS(total_nodes,m,n_near,x,weighting,lambda,n_loc,loc_index,kin,cor,x_c,region,ra,nearest_b,nearest_i)
      call MLS_zb(total_nodes,m,n_near,weighting,lambda,n_loc,loc_index,zb,zb_x,kin,nearest_i)

      step=0
      call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
      call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step) &
              ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))

      step=1  !在第０時間步和第１時間步之間插入 step=1/2，增加模式穩定
!      dt=dt/2
!      wavemaker_uh_bc = (wavemaker_uhArr(0) + wavemaker_uhArr(1)) / 2
!      wavemaker_h_bc = (wavemaker_hArr(0) + wavemaker_hArr(1)) / 2

      call predict(total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz,h_min,coriolis,uh_xx,u_max,kin)
!print*, 'after predict'
!print*, h(2494,step),h_x(2494,step),zb(2792)
      call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
      call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
      call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step) &
              ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))
!print*, 'after MLS'
!print*, h(2494,step),h_x(2494,step)
      do i_correct = 1,n_correct

          call calculate_q(total_nodes,step,dt,h,w_velocity,q)
          call MLS_q(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,q,q_x,h(:,step),h_min,zb,kin)
          call correct(i_correct,n_correct,total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz&
                        ,h_min,coriolis,uh_xx,u_max,kin,q,q_x)
!print*, 'after correct'
!print*, h(2494,step),h_x(2494,step)

          call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
          call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
          call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step)&
                  ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))


      end do


!      step=2  !其實是第１時間步，step=1
!      wavemaker_uh_bc = wavemaker_uhArr(1)
!      wavemaker_h_bc = wavemaker_hArr(1)
!      !write(*,'(a5,i5)') 'step=',step/2
!      call predict(total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz,h_min,coriolis,uh_xx,u_max,kin)
!!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
!      call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
!      call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
!      call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step) &
!              ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))
!!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
!      do i_correct = 1,n_correct
!
!          call calculate_q(total_nodes,step,dt,h,w_velocity,q)
!          call MLS_q(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,q,q_x,h(:,step),h_min,zb,kin)
!!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
!          call correct(i_correct,n_correct,total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz &
!                        ,h_min,coriolis,uh_xx,u_max,kin,q,q_x)
!!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
!          call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
!          call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
!          call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step) &
!                  ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))
!
!      end do

      ! for compare with experiment
      compareWithexperiment(1) = ( h(gaugePosition,step) + zb(gaugePosition) ) / h0

      do j=1,total_nodes

        h(j,step-1)=h(j,step)
        uh(j,step-1)=uh(j,step)
        h_x(j,step-1)=h_x(j,step)
        uh_x(j,step-1)=uh_x(j,step)
        u2h_x(j,step-1)=u2h_x(j,step)
        uh_xx(j,step-1)=uh_xx(j,step)
        w_velocity(j,step-1) = w_velocity(j,step)

        ! data for animation
        animation(j) = h(j,step) + zb(j)

      end do


      !dt=dt*2
      step=2 !計算第２時間步以後，即step>=2

      do step_n=2,total_steps  !新加一個計數index, step_n

        wavemaker_uh_bc = wavemaker_uhArr(step_n)
        wavemaker_h_bc = wavemaker_hArr(step_n)

        write(*,'(a5,i5)') 'step=',step_n   ! print the step
        call predict(total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz,h_min,coriolis,uh_xx,u_max,kin)
!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
        call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
        call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
        call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step)&
                ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))
!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
        do i_correct =1,n_correct

            call calculate_q(total_nodes,step,dt,h,w_velocity,q)
            call MLS_q(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,q,q_x,h(:,step),h_min,zb,kin)
!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
            call correct(i_correct,n_correct,total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz &
                            ,h_min,coriolis,uh_xx,u_max,kin,q,q_x)
!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
            call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
            call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
            call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step) &
                    ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))

        end do

        ! for compare with experiment
        compareWithexperiment(step_n) = ( h(gaugePosition,step) + zb(gaugePosition) ) / h0


2       do j=1,total_nodes
          h(j,0)=h(j,1)
          uh(j,0)=uh(j,1)
          h_x(j,0)=h_x(j,1)
          uh_x(j,0)=uh_x(j,1)
          u2h_x(j,0)=u2h_x(j,1)
          uh_xx(j,0)=uh_xx(j,1)
          w_velocity(j,0) = w_velocity(j,1)

          h(j,1)=h(j,2)
          uh(j,1)=uh(j,2)
          h_x(j,1)=h_x(j,2)
          uh_x(j,1)=uh_x(j,2)
          u2h_x(j,1)=u2h_x(j,2)
          uh_xx(j,1)=uh_xx(j,2)
          w_velocity(j,1) = w_velocity(j,2)

          ! data for animation
          animation((step_n - 1) * total_nodes + j) = h(j,step) + zb(j)    ! 1376 2902

        end do



        ! step data output

        if(step_n.eq.total_steps)then

          open(3,file='C:\Users\laserlab04\Desktop\zeta.txt')

          do j = 1,total_nodes
              write(3,'(11f16.5)')x(j) , h(j,step)+zb(j)
          end do

          close(3)


        end if


     end do

    ! output experiment data
!    open(3,file='C:\Users\laserlab04\Desktop\compareWithexperiment.txt')
!    do j = 1,total_steps
!        write(3,'(11f16.5)')compareWithexperiment(j)
!    end do
!    close(3)

    ! data for animation
!    open(3,file='C:\Users\laserlab04\Desktop\animation.txt')
!    do j = 1,total_nodes * total_steps
!
!        write(3,'(11f16.5)')animation(j)
!
!    end do
!    close(3)

    end

	! function ------------------------

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

    subroutine local_poly(x,p,m)
      implicit none
      integer m
      real*8 x,p(1:m)
      p(1)=1.0d0
      p(2)=x
      p(3)=x**2
    end

    subroutine KNN(N,x,d_near,n_near,cor,x_c,region,star)
      implicit none
      integer N,i,j,k,n_near,cor(1:N),region(1:N),star
      real*8 x(1:N),r,d_near(1:n_near),x_c(1:N)
  !  一開始最近距離都設1.d99
      do i=1,n_near
        d_near(i)=1.d99
      end do
  !  開始搜尋
      do j=1,N
        if(cor(star)*cor(j) .lt. 0) then
          if((x_c(star) .eq. x_c(j))) then
            if(((x_c(star)-x(star))-(x(j)-x(star)))*cor(star) .lt. 0) then
              goto 1  !有被遮蔽到的點要排除
            end if
          end if
        end if
        if(abs(region(star)-region(j)) .gt. 1) then
          goto 1  !不在同一區的點也要排除
        end if
        if(abs(x(j)-x(star)) .gt. 1.d5) goto 1
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

    subroutine MLS(N,x,m,n_near,w,lambda,n_loc,loc_index,h,uh,h_x,uh_x,u2h_x,zb,zb_x,h_min &
        ,uh_xx,u_max,ra,kin,w_velocity)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),wet_neighbors(1:N),wet(1:N),kin(1:N)&
             ,shore(1:N),n_nearest_wet
      real*8 x(1:N),y(1:N),h(1:N),uh(1:N),new_h(1:N),new_uh(1:N)&
            ,h_x(1:N),uh_x(1:N),u2h_x(1:N)&
            ,lambda(1:N,1:m,1:n_near*5),w(1:N,1:n_near*5),zb(1:N),zb_x(1:N),h_min &
            ,uh_xx(1:N),theta(1:N),eta(1:N),ww &
            ,smooth_eta(1:N),smooth_uh(1:N),err_eta(1:N),err_uh(1:N),err_vh(1:N)
      real*8 alpha(1:m),beta(1:n_near*5),p(1:m) &
            ,a(1:n_near*5,1:m),at(1:m,1:n_near*5),ata(1:m,1:m),ata_inv(1:m,1:m)
      real*8 x_loc,y_loc,r,u,v,u_max,abs_u,ra(1:N),r_nearest_wet,lambda_e(1:m,1:n_near*5)&
             ,w_e(1:n_near*5),r_nearest_dry(1:N)
      real*8 w_velocity(1:N),new_w_velocity(1:N),aa

      do j=1,N
        new_h(j)=h(j)
        new_uh(j)=uh(j)
        if(h(j) .le. h_min) then   !粗略判斷乾濕
          wet(j)=0
        else
          wet(j)=1
        end if
      end do
  !找出緊鄰濕點的乾點，並計算距最近乾點有多遠
      do j=1,N
        shore(j)=0  !先都設為0
        r_nearest_dry(j)=1.d99 !先都設為很大
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x(j)
          r=(x_loc**2)**.5d0
          if(r .lt. ra(j)*4.d0/3.d0) then                                      !緊鄰的點之中，
            if(wet(j) .eq. 0 .and. wet(loc_index(j,k)) .eq. 1) then            !有濕點，且這一點本身是乾點
              if(zb(j)+h(j) .lt. zb(loc_index(j,k))+h(loc_index(j,k))) then    !而本身乾點之高程比緊鄰濕點的水位還低
                shore(j)=1                                                     !則它是緊鄰濕點的乾點，要當成濕點處理
              end if
            end if
          end if
          if(r .lt. r_nearest_dry(j) .and. wet(loc_index(j,k)) .eq. 0) then    !緊鄰的點之中，是乾點且距離比到最近乾點還小
            r_nearest_dry(j)=r                                                 !則更新到最近乾點的距離
          end if
        end do
      end do
  !重新定義濕點
      do j=1,N
        wet(j)=wet(j)+shore(j) !本是濕點，則仍為1；若是乾點，但要當濕點處理，則0+1=1；乾點距濕點很遠，則為0
      end do
  !計算數量及比率
      do j=1,N
        wet_neighbors(j)=0
        do k=1,n_loc(j)
          if(wet(loc_index(j,k)) .eq. 1) then
            wet_neighbors(j)=wet_neighbors(j)+1
          end if
        end do
        theta(j)=1.d0*wet_neighbors(j)/n_loc(j) ! theta 代表周圍濕點的比率，若為1則表示周圍都是濕點，若為0則表示周圍都是乾點

      end do

  !給乾點定義eta值。乾點的eta並不存在，用最近濕點的eta值當它的eta。若最近濕點有兩個以上，則取平均。若範圍內完全無濕點，則用底床高程當eta值。
      do j=1,N
        if(theta(j) .eq. 0) then
          eta(j)=1.d99  !若範圍內完全無濕點，則設eta=1.99，在計算eta偏導數時，這類的點要排除，即要重新建構局部矩陣
          ! change this part to eta(j) = zb(j)
          !eta(j) = zb(j)
          goto 11
        end if
        if(wet(j) .eq. 1) then
          eta(j)=zb(j)+h(j)  !若本身是濕點，則用eta=zb+h
          goto 11
        end if
        eta(j)=0   !剩下的就是乾點
        ww=0
        do k=1,n_loc(j)
          if(wet(loc_index(j,k)) .eq. 1) then
            ww=ww+w(j,k)**2
            eta(j)=eta(j)+(h(loc_index(j,k))+zb(loc_index(j,k)))*w(j,k)**2
          end if
        end do
        eta(j)=eta(j)/ww
11    end do
  !開始計算
    ! h的部分
      do j=1,N
        if(kin(j) .eq. -1) goto 21            !ghost point用不到各物理量的偏導數。
        if(wet(j) .eq. 0 .or. theta(j) .lt. .15d0) then !若這一點是乾點，或範圍內濕點數量很少，則設為不流動。
          h_x(j)=-zb_x(j)
          goto 21
        end if
      !計算beta_eta
        do k=1,n_loc(j)
          if(eta(loc_index(j,k)) .gt. 1.d98) then
            goto 22 !若範圍內有eta=1.d99的點，則要重新建構局部矩陣
          end if
          !beta(k)=w(j,k)*h(loc_index(j,k))
          beta(k)=w(j,k)*eta(loc_index(j,k))
        end do
      !計算alpha_eta
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !輸出 h,h_x,h_y
        !new_h(j)=alpha(1)  ! h=eta-zb
        !new_h(j)=h(j)
        !h_x(j)=alpha(2)
         new_h(j)=eta(j) - zb(j)
         h_x(j)=alpha(2) - zb_x(j)

        goto 21
      !重新建構局部矩陣
      !計算A矩陣的內容
22      do k=1,n_loc(j)
          if(eta(loc_index(j,k)) .gt. 1.d98) then
            w_e(k)=0 !若eta=1.d99，則排除，設權重為0
          else
            w_e(k)=w(j,k) !不然就照原本的權重
          end if
          x_loc=x(loc_index(j,k))-x(j)
          call local_poly(x_loc,p,m) !計算Local Polynomial
          do i=1,m
            a(k,i)=w_e(k)*p(i)  ! 算eta時，只拿eta不等於1.d99的點來算。
          end do
        end do
      !產生A矩陣的轉置矩陣
        do k=1,n_loc(j)
          do i=1,m
            at(i,k)=a(k,i)
          end do
        end do
      !產生Psuudo矩陣
        do i=1,m
          do l=1,m
            ata(i,l)=0
            do k=1,n_loc(j)
              ata(i,l)=ata(i,l)+at(i,k)*a(k,l)    !atxa=at*a
            end do
          end do
        end do
        call matrix_inv(m,ata,ata_inv)
      !產生lambda矩陣
        do i=1,m
          do l=1,n_loc(j)
            lambda_e(i,l)=0
            do k=1,m
              lambda_e(i,l)=lambda_e(i,l)+ata_inv(i,k)*at(k,l)
            end do
          end do
        end do
      !計算beta_eta
        do k=1,n_loc(j)
          beta(k)=w_e(k)*h(loc_index(j,k))
        end do
      !計算alpha_eta
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda_e(i,k)*beta(k)
          end do
        end do
      !輸出 h,h_x,h_y
        !new_h(j)=alpha(1) ! h=eta-zb
        new_h(j)=h(j)
        h_x(j)=alpha(2)

21    end do


    ! uh,u2h
      do j=1,N
        if(kin(j) .eq. -1) goto 31            !ghost point用不到各物理量的偏導數。
        if(wet(j) .eq. 0 .or. theta(j) .lt. .15d0) then !若這一點是乾點，或範圍內濕點數量很少，則設為不流動。
          new_uh(j)=0
          uh_x(j)=0
          uh_xx(j)=0
          u2h_x(j)=0
          goto 31
        end if
      !計算beta_uh
        do k=1,n_loc(j)
          beta(k)=w(j,k)*uh(loc_index(j,k))
        end do
      !計算alpha_uh
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !輸出 uh,uh_x,uh_y
        !new_uh(j)=alpha(1)
        new_uh(j)=uh(j)
        uh_x(j)=alpha(2)
        uh_xx(j)=2*alpha(4)
      !計算beta_u2h
        do k=1,n_loc(j)
          if(h(loc_index(j,k)) .le. h_min) then
            beta(k)=0
          else
            beta(k)=w(j,k)*uh(loc_index(j,k))**2/h(loc_index(j,k))
          end if
        end do
      !計算alpha_u2h
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !輸出 u2h_x
        u2h_x(j)=alpha(2)
31    end do

!-----------------------------------------------------------------------------
    !檢查是否有水深為負值或流速過大
      ww=3.5d0
      do j=1,N
        if(kin(j) .eq. -1) goto 32

        if(r_nearest_dry(j) .gt. ww*ra(j)) then !到最近乾點距離超過ww倍的ra
          h(j)=new_h(j)
          uh(j)=new_uh(j)
        else
          h(j)=new_h(j)*(r_nearest_dry(j)/(ww*ra(j)))+h(j)*(1-r_nearest_dry(j)/(ww*ra(j)))
          uh(j)=new_uh(j)
        end if

        if(h(j) .lt. 0) then
          h(j)=0
        end if

        if(h(j) .le. h_min) then
          uh(j)=0
          goto 32
        end if

        u=uh(j)/h(j)
        abs_u=(u**2)**.5d0

!        if(abs_u .gt. u_max) then
!          uh(j)=uh(j)*u_max/abs_u
!        end if

!------------------------------------------------------------------

32    end do

!      print*, 'MLS '
!      print*, h(2792),h_x(2792),zb(2792),zb_x(2792)
      !print*, 'MLS uh uh_x'
      !print*, uh(2792),uh_x(2792)

        do j = 1,N

            if(wet(j) .eq. 0 .or. theta(j) .lt. .15d0 ) then

               w_velocity(j) = 0

            else

               w_velocity(j) = 0.5 * ( -uh_x(j) + ( uh(j) / h(j) ) * (h_x(j) + 2 * zb_x(j) ) )

            end if

        end do

!      ! smooth w
!           do i=1,8
!            do j=1,N
!              aa=0
!              new_w_velocity(j)=0
!              do k=1,n_loc(j)
!                aa=aa+w(j,k)**2
!                new_w_velocity(j)=new_w_velocity(j)+w(j,k)**2*w_velocity(loc_index(j,k))
!              end do
!              new_w_velocity(j)=new_w_velocity(j)/aa
!            end do
!            do j=1,N
!              w_velocity(j)=new_w_velocity(j)
!            end do
!          end do

!      !開始計算
!          do j=1,N
!          !計算beta_w
!            do k=1,n_loc(j)
!              beta(k)=w(j,k)*w_velocity(loc_index(j,k))
!            end do
!          !計算alpha_w
!            do i=1,m
!              alpha(i)=0.d0
!              do k=1,n_loc(j)
!                alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
!              end do
!            end do
!          !輸出 smooth_w
!            new_w_velocity(j)=alpha(1)
!            new_w_velocity(j)=w_velocity(j)
!          end do
!      !代換
!          do j=1,N
!            w_velocity(j)=new_w_velocity(j)
!          end do

    end


    subroutine predict(total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz0,h_min,coriolis,uh_xx,u_max,kin)
      implicit none
      integer j,total_nodes,step,kin(1:total_nodes)
      real*8 h(1:total_nodes,0:2),uh(1:total_nodes,0:2),h_x(1:total_nodes,0:2),uh_x(1:total_nodes,0:2) &
            ,u2h_x(1:total_nodes,0:2),dt,g,zb_x(1:total_nodes),Cz,Cz0,u,abs_u,h_min,coriolis(1:total_nodes) &
            ,Eddy,uh_xx(1:total_nodes,0:2),u_max

      Eddy=0.d0

      do j=1,total_nodes
        if(kin(j) .eq. -1) goto 1

        if(step .eq. 1) then

          ! h------
          h(j,step)=h(j,step-1)-dt * uh_x(j,step-1)

          if(h(j,step) .lt. 0) then
            h(j,step)=0
          end if

          ! uh------
          if(h(j,step) .le. h_min) then
            uh(j,step)=0
            goto 1
          end if

          ! step - 1
          if(h(j,step-1) .le. h_min) then
            u=0
          else
            u=uh(j,step-1)/h(j,step-1)
          end if

          uh(j,step)=uh(j,step-1)-dt * u2h_x(j,step-1)                             !convective terms
          uh(j,step)=uh(j,step)-dt*g*h(j,step-1)*(zb_x(j)+h_x(j,step-1))           !variable bottom elevation and gravity terms

        else

          ! h------
          h(j,step)=h(j,step-1)-dt * (1.5d0 * uh_x(j,step-1) - .5d0 * uh_x(j,step-2))

          if(h(j,step) .lt. 0) then
            h(j,step)=0
          end if

          ! uh-----
          if(h(j,step) .le. h_min) then
            uh(j,step)=0
            goto 1
          end if

          ! step - 1
          if(h(j,step-1) .le. h_min) then
            u=0
          else
            u=1.5d0*uh(j,step-1)/h(j,step-1)
          end if

          ! step - 2
          if(h(j,step-2) .le. h_min) then
            u=u
          else
            u=u-.5d0*uh(j,step-2)/h(j,step-2)
          end if


          uh(j,step)=uh(j,step-1)-dt*( 1.5d0 * u2h_x(j,step-1) - .5d0 * u2h_x(j,step-2) )
          uh(j,step)=uh(j,step)-dt*g*(1.5d0 * h(j,step-1) * (zb_x(j)+h_x(j,step-1)) &
                                    - .5d0 * h(j,step-2) * (zb_x(j)+h_x(j,step-2)))

        end if

        u=uh(j,step)/h(j,step)


1     end do

!    print*, 'in predict'
!    print*, zb_x(2792),h_x(2792,step-1)
!    print*, zb_x(2792),h_x(2792,step-2)

    end


    subroutine correct(i_correct,n_correct,total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz0,h_min&
        ,coriolis,uh_xx,u_max,kin,q,q_x)
      implicit none
      integer j,total_nodes,step,kin(1:total_nodes),i_correct,n_correct
      real*8 h(1:total_nodes,0:2),uh(1:total_nodes,0:2),h_x(1:total_nodes,0:2),uh_x(1:total_nodes,0:2)&
            ,u2h_x(1:total_nodes,0:2),dt,g,zb_x(1:total_nodes),Cz,Cz0,u,v,abs_u,h_min,coriolis(1:total_nodes) &
            ,Eddy,uh_xx(1:total_nodes,0:2),u_max
      real*8 q(1:total_nodes),q_x(1:total_nodes),rr

      Eddy=0.d0
      rr=1.d0*i_correct/n_correct

      do j=1,total_nodes
        if(kin(j) .eq. -1) goto 1
        ! h------
        h(j,step)=h(j,step-1)-dt * 0.5 * (uh_x(j,step)+uh_x(j,step-1))

        if(h(j,step) .lt. 0) then
          h(j,step)=0
        end if

        ! uh------
        if(h(j,step) .le. h_min) then
          uh(j,step)=0
          goto 1
        end if
        u=uh(j,step)/h(j,step)

        if(h(j,step-1) .gt. h_min) then
          u=u+uh(j,step-1)/h(j,step-1)
        end if

        u=u/2

        uh(j,step)=uh(j,step-1)- 0.5 * dt*(u2h_x(j,step)+u2h_x(j,step-1))                       !convective terms

        uh(j,step)=uh(j,step)- 0.5 * dt*g*(h(j,step)*(zb_x(j)+h_x(j,step)) &
                                   +h(j,step-1)*(zb_x(j)+h_x(j,step-1)))                 !variable bottom elevation and gravity terms

        ! non-hydrostatic terms


        uh(j,step) = uh(j,step) - 0.5 * dt * 0.5 * (h(j,step-1) + h(j,step)) * q_x(j)

        uh(j,step) = uh(j,step) - 0.5 * dt * q(j) * 0.5 * ( h_x(j,step-1) + h_x(j,step) + 2*zb_x(j) )


        u=uh(j,step)/h(j,step)


1     end do


!        print*, 'in correct'
!        print*, zb_x(2792),h_x(2792,step-1)
!        print*, zb_x(2792),h_x(2792,step)


    end


    subroutine NF_bc(N,h,uh,kin,wavemaker_uh_bc,wavemaker_h_bc)
      implicit none
      integer j,N,kin(1:N)
      real*8 h(1:N),uh(1:N),wavemaker_uh_bc,wavemaker_h_bc

      do j=1,N
        ! no flow bc

        if(kin(j) .eq. 1) then
            uh(j) = 0
            h(j) = 0
        end if

        ! wave maker bc

        if(kin(j) .eq. 2) then
            uh(j) = 0
            h(j) = 0.218
            !uh(j) = wavemaker_uh_bc
            !h(j) = wavemaker_h_bc
        end if

      end do
    end


    subroutine NS_bc(N,uh,vh,kin)  ! not using this function
      implicit none
      integer j,N,kin(N)
      real*8 uh(1:N),vh(1:N)
      do j=1,N
        if(kin(j) .eq. 2) then
          uh(j)=0
          vh(j)=0
        end if
      end do
    end


!    subroutine FF_bc(N,n_x,h,uh,kin,h_min)  ! not using this function
!      implicit none
!      integer j,N,kin(1:N)
!      real*8 n_x(1:N),h(1:N),uh(1:N),u,v,abs_u,cr_u,h_min
!      do j=1,N
!        if(kin(j) .eq. 3) then
!          if(h(j) .gt. h_min) then
!            u=uh(j)/h(j)
!            v=vh(j)/h(j)  !----v can be deleted----
!            abs_u=(u**2+v**2)**.5  !----abs_u = abs(u)----
!            cr_u=(9.81*h(j))**.5
!            if(abs_u .lt. cr_u) then
!              uh(j)=n_x(j)*cr_u*h(j)
!              vh(j)=n_y(j)*cr_u*h(j)
!            end if
!          end if
!        end if
!      end do
!
!      do j = 1,N
!
!        if(kin(j) .eq. 2)then
!
!        end if
!
!      end do
!
!    end


    subroutine eta_bc(N,h,zb,kin,eta,n_eta_BC)
      implicit none
      integer j,N,kin(N),n_eta_BC
      real*8 h(1:N),zb(1:N),eta(1001:1000+n_eta_BC)
      do j=1,N
        if((kin(j) .gt. 1.d3) .and. (kin(j) .le. 5.d3)) then
          h(j)=eta(kin(j))-zb(j)
        end if
      end do
    end


    subroutine flux_bc(N,uh,vh,kin,flux_x,flux_y,n_flux_BC)  ! not using this function
      implicit none
      integer j,N,kin(N),n_flux_BC
      real*8 uh(1:N),vh(1:N),flux_x(5001:5000+n_flux_BC),flux_y(5001:5000+n_flux_BC)
      do j=1,N
        if(kin(j) .gt. 5.d3) then
          uh(j)=flux_x(kin(j))
          vh(j)=flux_y(kin(j)) !----vh----
        end if
      end do
    end


    subroutine MLS_zb(N,m,n_near,w,lambda,n_loc,loc_index,zb,zb_x,kin,nearest_i)
      implicit none
      integer i,j,k,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),kin(1:N),nearest_i(1:N)
      real*8 zb(1:N),zb_x(1:N),new_zb(1:N),lambda(1:N,1:m,1:n_near*5),w(1:N,1:n_near*5)
      real*8 alpha(1:m),beta(1:n_near*5)
! 先把ghost point的zb換成對應內部點之zb
      do j=1,N
        if(kin(j) .eq. -1) then
          zb(j)=zb(nearest_i(j))
        end if
      end do
!開始計算
      do j=1,N
  !計算beta_zb
        do k=1,n_loc(j)
          beta(k)=w(j,k)*zb(loc_index(j,k))
        end do
  !計算alpha_zb
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
  !輸出 zb,zb_x,zb_y
        !new_zb(j)=alpha(1)
        new_zb(j)=zb(j)

        if(j .eq. 3000)then
          zb_x(j)=0
        else
          zb_x(j)=alpha(2)
        end if

      end do
      do j=1,N
        zb(j)=new_zb(j)
      end do
    end


    subroutine pre_MLS(N,m,n_near,x,w,lambda,n_loc,loc_index,kin,cor,x_c,region,ra,nearest_b,nearest_i)
      implicit none
      integer i,j,k,l,m,N,n_near,kin(1:N),n_loc(1:N),loc_index(1:N,1:n_near*5),cor(1:N),region(1:N) &
             ,nearest_b(1:N),nearest_i(1:N)
      real*8 x(1:N),w(1:N,1:n_near*5),lambda(1:N,1:m,1:n_near*5),d_near(1:n_near),p(1:m),x_c(1:N)
      real*8 r,x_loc,ro,epsilon,ra(1:N),r_nearest_b
      real*8 a(1:n_near*5,1:m),at(1:m,1:n_near*5),ata(1:m,1:m),ata_inv(1:m,1:m),x_i
      epsilon=-10d0
  !開始做Collocation
      do j=1,N
  !找第n_near最近點的距離
        call KNN(N,x,d_near,n_near,cor,x_c,region,j)
        ra(j)=(d_near(2)+d_near(3)+d_near(4)+d_near(5))/4.d0 ! ra 大概是local nodal spacing
        ro=max(ra(j)*4.5d0,d_near(n_near))*1.001
  !搜點並計算權重
        n_loc(j)=0
        do l=1,N !逐一搜點
          if(cor(l)*cor(j) .lt. 0) then
            if((x_c(l) .eq. x_c(j))) then
              if(((x_c(j)-x(j))-(x(l)-x(j)))*cor(j) .lt. 0) then   ! might cause problems
                goto 1 !有被遮蔽到的點要排除
              end if
            end if
          end if
          if(abs(region(l)-region(j)) .gt. 1) then
            goto 1 !不在同一區的點也要排除
          end if
          x_loc=x(l)-x(j) !計算相對坐標
          r=abs(x_loc) !計算到focused node之間的距離  !----r = x_loc----
          if(r .ge. ro) goto 1  !若距focused node太遠，也要排除
          n_loc(j)=n_loc(j)+1
          loc_index(j,n_loc(j))=l !這裡的local_index就是論文裡的k。即第l個node在第j個局部範圍內，其編號為k
          w(j,n_loc(j))=(1-r/ro)**3.65d0 !abs((exp(epsilon*(r/ro)**2)-exp(epsilon))/(1.d0-exp(epsilon)))**.5d0 !計算權重係數
1       end do
      end do
! 搜尋ghost point所對應的邊界點
      do j=1,N
        if(kin(j) .ge. 0) goto 2
        r_nearest_b=1.d99
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x(j) !計算相對坐標
          r=abs(x_loc)  !計算到focused node之間的距離  !----r = x_loc----
          if(kin(loc_index(j,k)) .gt. 0) then !局部第k點是邊界點，內部點kin是0，ghost point的kin是-1，不是內部也不是ghost就是邊界點
            if(r .lt. r_nearest_b) then !距離比原本的r_nearest_b還小
              r_nearest_b=r
              nearest_b(j)=loc_index(j,k) !得到最近邊界點的局部編號
            end if
          end if
        end do
2     end do
! 搜尋ghost point所對應的內部點
      do j=1,N
        if(kin(j) .ge. 0) goto 3
        x_i=2*x(nearest_b(j))-x(j) !x_i, y_i為對應內部點之位置
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x_i  !先計算
          r=abs(x_loc)
          if(r .lt. 1.d-6) then
            nearest_i(j)=loc_index(j,k) !若r近乎0，則它就是ghost point所對應之內部點
            goto 3
          end if
        end do
3     end do
!建構局部矩陣
      do j=1,N
  !計算A矩陣的內容
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x(j)  !計算
          call local_poly(x_loc,p,m) !計算Local Polynomial
          do i=1,m
            a(k,i)=w(j,k)*p(i)
          end do
        end do
  !產生A矩陣的轉置矩陣
        do k=1,n_loc(j)
          do i=1,m
            at(i,k)=a(k,i)
          end do
        end do
  !產生Psuudo矩陣
        do i=1,m
          do l=1,m
            ata(i,l)=0
            do k=1,n_loc(j)
              ata(i,l)=ata(i,l)+at(i,k)*a(k,l)    !atxa=at*a
            end do
          end do
        end do
  !產生Psuudo矩陣的反矩陣
        call matrix_inv(m,ata,ata_inv)
  !產生lambda矩陣
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


    subroutine og_bc(N,n_near,n_loc,loc_index,kin,x,y,h,uh,vh,ra)  ! not using this function
      implicit none
      integer j,k,N,kin(N),n_loc(1:N),loc_index(1:N,1:n_near*5),n_near,n_count
      real*8 h(1:N),uh(1:N),vh(1:N),x(1:N),y(1:N),ra(1:N),r,x_loc,y_loc
      do j=1,N
        if(kin(j) .eq. 4) then
          h(j)=0
          uh(j)=0
          vh(j)=0           !----vh----
          n_count=0
          do k=1,n_loc(j)
            x_loc=x(loc_index(j,k))-x(j)
            y_loc=y(loc_index(j,k))-y(j)  !----y_loc----
            r=(x_loc**2+y_loc**2)**.5     !----r = x_loc----
            if(r .le. 1.5d0*ra(j)) then !緊臨
              if(kin(loc_index(j,k)) .ne. 4) then !非out going邊界
                n_count=n_count+1
                h(j)=h(j)+h(loc_index(j,k))
                uh(j)=uh(j)+uh(loc_index(j,k))
                vh(j)=vh(j)+vh(loc_index(j,k))  !----vh----
              end if
            end if
          end do
          h(j)=h(j)/n_count  !out going邊界點的h值等於緊臨非out going邊界點h值之平均
          uh(j)=uh(j)/n_count !out going邊界點的uh值等於緊臨非out going邊界點uh值之平均
          vh(j)=vh(j)/n_count !out going邊界點的vh值等於緊臨非out going邊界點vh值之平均 !----vh----
        end if
      end do
    end



    subroutine ghost_BC(N,kin,h,uh,nearest_b,nearest_i,x)
      implicit none
      integer j,N,kin(N),nearest_b(1:N),nearest_i(1:N)
      real*8 h(1:N),uh(1:N),x(1:N),n_x,n_abs,A
      do j=1,N
        if(kin(j) .eq. -1) then

          h(j)=h(nearest_i(j))
          n_x=x(j)-x(nearest_b(j))
          n_abs=(n_x**2)**.5d0
          n_x=n_x/n_abs
          A = -n_x**2
          uh(j) = A * uh(nearest_i(j))

        end if
      end do
    end


    subroutine calculate_q(total_nodes,step,dt,h,w_velocity,q)
      implicit none
      integer j,total_nodes,step
      real*8 h(1:total_nodes,0:2),w_velocity(1:total_nodes,0:2),q(1:total_nodes),dt,W_t
      do j=1,total_nodes
	    W_t = (w_velocity(j,step)-w_velocity(j,step-1)) / dt

        q(j) = 0.5 * (h(j,step)+h(j,step-1)) * W_t

      end do

    end

    subroutine MLS_q(N,x,m,n_near,w,lambda,n_loc,loc_index,q,q_x, h,h_min,zb,kin)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),n_smooth,kin(1:N)
      real*8 x(1:N),q(1:N),q_x(1:N),lambda(1:N,1:m,1:n_near*5),w(1:N,1:n_near*5)
      real*8 alpha(1:m),beta(1:n_near*5),smooth_q(1:N),aa
      ! part 1
      real*8 wet(1:N),h_min,h(1:N)
      ! part 2
      real*8 shore(1:N),r_nearest_dry(1:N),x_loc,zb(1:N), r, ra(1:N)
      ! part 3
      real*8 wet_neighbors(1:N),theta(1:N),eta(1:N),ww
      real*8 lambda_e(1:m,1:n_near*5),w_e(1:n_near*5)
      real*8 p(1:m),a(1:n_near*5,1:m),at(1:m,1:n_near*5),ata(1:m,1:m),ata_inv(1:m,1:m)

	  n_smooth=8

	! part 1
      do j=1,N
        if(h(j) .le. h_min) then   !粗略判斷乾濕
          wet(j)=0
        else
          wet(j)=1
        end if
      end do

    ! part 2
      !找出緊鄰濕點的乾點，並計算距最近乾點有多遠
      do j=1,N
        shore(j)=0  !先都設為0
        r_nearest_dry(j)=1.d99 !先都設為很大
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x(j)
          r=(x_loc**2)**.5d0
          if(r .lt. ra(j)*4.d0/3.d0) then                                      !緊鄰的點之中，
            if(wet(j) .eq. 0 .and. wet(loc_index(j,k)) .eq. 1) then            !有濕點，且這一點本身是乾點
              if(zb(j)+h(j) .lt. zb(loc_index(j,k))+h(loc_index(j,k))) then    !而本身乾點之高程比緊鄰濕點的水位還低
                shore(j)=1                                                     !則它是緊鄰濕點的乾點，要當成濕點處理
              end if
            end if
          end if
          if(r .lt. r_nearest_dry(j) .and. wet(loc_index(j,k)) .eq. 0) then    !緊鄰的點之中，是乾點且距離比到最近乾點還小
            r_nearest_dry(j)=r                                                 !則更新到最近乾點的距離
          end if
        end do
      end do

    !重新定義濕點
      do j=1,N
        wet(j)=wet(j)+shore(j) !本是濕點，則仍為1；若是乾點，但要當濕點處理，則0+1=1；乾點距濕點很遠，則為0
      end do

    ! part 3
      !計算數量及比率
      do j=1,N
        wet_neighbors(j)=0
        do k=1,n_loc(j)
          if(wet(loc_index(j,k)) .eq. 1) then
            wet_neighbors(j)=wet_neighbors(j)+1
          end if
        end do
        theta(j)=1.d0*wet_neighbors(j)/n_loc(j) ! theta 代表周圍濕點的比率，若為1則表示周圍都是濕點，若為0則表示周圍都是乾點
      end do
      !給乾點定義eta值。乾點的eta並不存在，用最近濕點的eta值當它的eta。若最近濕點有兩個以上，則取平均。若範圍內完全無濕點，則用底床高程當eta值。
      do j=1,N
        if(theta(j) .eq. 0) then
          eta(j)=1.d99  !若範圍內完全無濕點，則設eta=1.99，在計算eta偏導數時，這類的點要排除，即要重新建構局部矩陣
          goto 11
        end if
        if(wet(j) .eq. 1) then
          eta(j)=zb(j)+h(j)  !若本身是濕點，則用eta=zb+h
          goto 11
        end if
        eta(j)=0   !剩下的就是乾點
        ww=0
        do k=1,n_loc(j)
          if(wet(loc_index(j,k)) .eq. 1) then
            ww=ww+w(j,k)**2
            eta(j)=eta(j)+(h(loc_index(j,k))+zb(loc_index(j,k)))*w(j,k)**2
          end if
        end do
        eta(j)=eta(j)/ww
11    end do


    ! Q
    !先把Q平滑一下
    do i=1,n_smooth
      do j=1,N
	    aa=0
		smooth_q(j)=0
        do k=1,n_loc(j)
		  aa=aa+w(j,k)**2
          smooth_q(j)=smooth_q(j)+w(j,k)**2*q(loc_index(j,k))
        end do
        smooth_q(j)=smooth_q(j)/aa
	  end do
      do j=1,N
        q(j)=smooth_q(j)
	  end do
	end do

    ! Q的部分
      do j=1,N
        if(kin(j) .eq. -1) goto 21            !ghost point用不到各物理量的偏導數。
        if(wet(j) .eq. 0 .or. theta(j) .lt. .15d0) then !若這一點是乾點，或範圍內濕點數量很少，則設為不流動。
          q(j) = 0
          q_x(j) = 0
          goto 21
        end if
      !計算beta_q
        do k=1,n_loc(j)
          if(eta(loc_index(j,k)) .gt. 1.d98) then
            goto 22 !若範圍內有eta=1.d99的點，則要重新建構局部矩陣
          end if
          beta(k)=w(j,k)*q(loc_index(j,k))
        end do
      !計算alpha_q
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !輸出 q,q_x
        !q(j)=alpha(1)  ! h=eta-zb
        q_x(j)=alpha(2)
        goto 21
      !重新建構局部矩陣
      !計算A矩陣的內容
22      do k=1,n_loc(j)
          if(eta(loc_index(j,k)) .gt. 1.d98) then
            w_e(k)=0 !若eta=1.d99，則排除，設權重為0
          else
            w_e(k)=w(j,k) !不然就照原本的權重
          end if
          x_loc=x(loc_index(j,k))-x(j)
          call local_poly(x_loc,p,m) !計算Local Polynomial
          do i=1,m
            a(k,i)=w_e(k)*p(i)  ! 算eta時，只拿eta不等於1.d99的點來算。
          end do
        end do
      !產生A矩陣的轉置矩陣
        do k=1,n_loc(j)
          do i=1,m
            at(i,k)=a(k,i)
          end do
        end do
      !產生Psuudo矩陣
        do i=1,m
          do l=1,m
            ata(i,l)=0
            do k=1,n_loc(j)
              ata(i,l)=ata(i,l)+at(i,k)*a(k,l)    !atxa=at*a
            end do
          end do
        end do
        call matrix_inv(m,ata,ata_inv)
      !產生lambda矩陣
        do i=1,m
          do l=1,n_loc(j)
            lambda_e(i,l)=0
            do k=1,m
              lambda_e(i,l)=lambda_e(i,l)+ata_inv(i,k)*at(k,l)
            end do
          end do
        end do
      !計算beta_eta
        do k=1,n_loc(j)
          beta(k)=w_e(k)*q(loc_index(j,k))
        end do
      !計算alpha_eta
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda_e(i,k)*beta(k)
          end do
        end do
      !輸出 h,h_x,h_y
        !q(j)=alpha(1) ! h=eta-zb
        q_x(j)=alpha(2)
21    end do

    end

