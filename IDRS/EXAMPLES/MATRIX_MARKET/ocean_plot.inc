! Make plot?
   plot = get_parameter('-plot', default_plot ) 
   plot_solution = ( plot == 'solution' .or. plot == 'all' )
   
! Make gifs?
   gif = get_parameter('-gif', default_gif )
   if ( gif ) then                
      plot = get_parameter('-gif', default_plot ) 
      plot_solution = ( plot == 'solution' .or. plot == 'all' )
   end if
   
! Make jpeg?
   jpeg = get_parameter('-jpeg', default_jpeg )
   if ( jpeg ) then               
      plot = get_parameter('-jpeg', default_plot )
      plot_solution = ( plot == 'solution' .or. plot == 'all' )
   end if

   if ( plot_solution ) then
! Read matrix that maps to earth coordinates:
      write(map_file,'(a,a,i1,a)') trim(path),'STOMMEL/stommel', grid,'_map.mtx'

! Read the map matrix:
      Q = mm_matrix( trim(map_file), 'coo' )
      ngrid = Q%nrows                ! Number of local gridpoints
      allocate(solution(ngrid,nrhs))
      do irhs = 1, nrhs
         solution(:,irhs) = Q*x(:,irhs)
      end do

      nx = 360/grid+1
      ny = 180/grid
      if ( my_proc == 1 ) then
         open(30,file='solution.dat')
      else
         sync images (my_proc-1)
         open(30,file='solution.dat',position = 'append',Status='old')
      end if
      nb = (my_proc-1)*ngrid         ! Last equation on previous image
      ne = min(my_proc*ngrid,nx*ny)  ! Last equation on this image
      do igrid = nb+1, ne
         i = mod(igrid-1,nx)+1
         j = (igrid-1)/nx+1
         write(30,'(i3,1x,i3,*(1x,e9.2))') i, j, real(solution(igrid-nb,:),kind=rp)
         if ( mod(igrid,nx) == 0 ) write(30,*)
      end do
      close(30)
      if ( my_proc < num_procs ) sync images (my_proc+1)

      if ( my_proc == 1 ) then
         open(40,file='solution.plt')
         if ( gif )  then 
            write(40,'(a)') 'set terminal gif animate delay 1'
            write(40,'(a)') 'set output "solution.gif"'
         end if
         if ( jpeg ) then
            write(40,'(a)') 'set terminal jpeg'
            write(40,'(a)') 'set output "solution.jpeg"'
            write(40,'(a)') 'set multiplot layout 3,4'
         end if
         write(40,'(a)') 'unset key'
         write(40,'(a)') 'unset border'
         write(40,'(a)') 'unset colorbox'
         write(40,'(a)') 'unset tics'
         write(40,'(a,f5.2)') 'set size ratio ', 0.5
         write(40,'(a)') 'set view projection xy'
         write(40,'(a)') 'set contourfill auto 20'
         do irhs = 1, nrhs
            write(40,'(a,a,a)') 'set title "', trim(month(irhs)), '"'
            write(40,'(a,i2,a)') 'splot "solution.dat" using 1:2:',irhs+2, ' with contourfill'
            if ( .not. gif .and. .not. jpeg ) write(40,'(a)') 'pause 1'
         end do
         close(40)
         call execute_command_line('gnuplot solution.plt')
      end if
   end if
 
