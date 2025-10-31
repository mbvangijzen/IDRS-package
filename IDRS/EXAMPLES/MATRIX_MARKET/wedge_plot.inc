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
      ngrid = nrows                      ! Number of local gridpoints
      nx = 2**(grid-1)*6+1
      ny = 2**(grid-1)*10+1
      if ( my_proc == 1 ) then
         open(30,file='solution.dat')
      else
         sync images (my_proc-1)
         open(30,file='solution.dat',position = 'append',Status='old')
      end if
      nb = (my_proc-1)*ngrid           ! Last gridpoint on previous image
      ne = min(my_proc*ngrid,nx*ny)    ! Last gridpoint on this image
      do igrid = nb+1, ne
         i = mod(igrid-1,nx)+1
         j = (igrid-1)/nx+1
         write(30,'(i3,1x,i3,*(1x,e12.5))') i, ny-j+1, real(x(igrid-nb,:),kind=rp)
         if ( mod(igrid,nx) == 0 ) write(30,*)
      end do
      close(30)
      if ( my_proc < num_procs ) sync images (my_proc+1)

      if ( my_proc == 1 ) then
         open(40,file='solution.plt')  
         if ( gif )  then
            write(40,'(a)') 'set terminal gif'
            write(40,'(a)') 'set output "solution.gif"'
         elseif ( jpeg ) then
            write(40,'(a)') 'set terminal jpeg'
            write(40,'(a)') 'set output "solution.jpeg"'
         end if
         write(40,'(a)') 'unset key'    
         write(40,'(a)') 'unset border'
         write(40,'(a)') 'unset colorbox'
         write(40,'(a)') 'unset tics'
         write(40,'(a)') 'set palet gray'
         select case (grid)
            case (1:3)
               i = 1
               j = grid
            case (4)
               i = 2
               j = 2
            case (5)
               i = 1
               j = 5
            case (6)
               i = 2
               j = 3
         end select
         write(40,'(a,i1,a,i1)') 'set multiplot layout ',i,',',j
         write(40,'(a,f5.2)') 'set size ratio ', 10./6.
         do isigma = 1, Nsigma
            write(40,'(a,f5.2,a)') 'set title "Frequency = ', freq(isigma), ' Hz"'
            write(40,'(a,i1,a)') 'plot "solution.dat" using 1:2:',isigma+2, ' with image'
         end do
         close(40)
         call execute_command_line('gnuplot -p solution.plt')
      end if
   end if
