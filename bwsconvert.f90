!*******************************************************************n*

program bws_convert_coords

!********************************************************************
!
! Originally written by MDE to convert coordinates to the AVS
! format. Simple bits modified by MW and MJC.
!
! Stuttgart (18.5.95): Modified to write a file in the `bws' format.
!
! 25.03.10 converted to fortran 90 by BJM
!
! 27.06.12 BJM: 
!       Dipole plotting removed. 
!       Rewritten to use user-defined types
!       Changed to write out atom flag status for passing additional
!       information to xbs (also updated)
!
!********************************************************************

implicit none

type atoms
    double precision, dimension(3) :: r
    integer :: ntype
    integer :: flag
end type atoms

type pairs
    character(len=5) :: name1
    character(len=5) :: name2
    double precision :: minl
    double precision :: maxl
    double precision :: bondgrey
    double precision :: bondradius
end type pairs

type species
    integer :: num
    type(atoms), dimension(:), allocatable :: atom
    double precision :: radius
    double precision, dimension(3) :: rgb
    character(len=5) :: name
    character(len=25) :: flag_filename
    logical :: read_flags
    integer :: fileflags
end type species

double precision :: h(3,3)
double precision :: boxlen(3)
double precision, dimension(3) :: rtemp

integer :: nframes, nspec, npairs
integer :: i, j, k
integer :: nskip
integer :: filein, fileout, filebox
integer :: err

character(len=80) :: filename, filename_out, filename_box, mvin, mvout

logical :: write_movie

type(species), allocatable, dimension(:) :: spec
type(pairs), allocatable, dimension(:) :: pair

interface
    
    function cell_to_lab( r, h )
        double precision, dimension(3), intent(in) :: r
        double precision, dimension(3,3), intent(in) :: h
        double precision, dimension(3) :: cell_to_lab
    end function cell_to_lab

end interface 

! Read in control variables.
open (file='bwsconvert.inpt', status='old', newunit=filein)

read (filein,*) nspec
allocate(spec(nspec), stat=err)
if (err /= 0) print *, "memory allocation failed for species"

do i=1, nspec
    associate( this_spec => spec(i) )
        read(filein,*) this_spec%num
        allocate( this_spec%atom( this_spec%num ), stat=err )
        if (err /= 0) print *, "memory allocation failed for atoms"
        read(filein,*) this_spec%radius
        read(filein,*) this_spec%rgb
        read(filein,'(a)') this_spec%name
        read(filein,*) this_spec%read_flags
        if (this_spec%read_flags) read(filein,'(a)') this_spec%flag_filename
    end associate
end do

read(filein,*) npairs
allocate( pair(npairs) )
if (err /= 0) print *, "memory allocation failed for pairs"
print *,'npairs',npairs
do i=1, npairs
    read( filein, '(a)' ) pair(i)%name1
    read( filein, '(a)' ) pair(i)%name2
    print *,'name1(',i,') = ', pair(i)%name1
    print *,'name2(',i,') = ', pair(i)%name2
    read( filein, *) pair(i)%minl
    read( filein, *) pair(i)%maxl
    read( filein, *) pair(i)%bondgrey
    read( filein, *) pair(i)%bondradius
end do

read( filein, '(a)') filename
read( filein, '(a)') filename_box
read( filein, '(a)') filename_out
read( filein,* ) write_movie

if (write_movie) then
    read( filein, *) nskip
    read( filein, *) nframes
    read( filein, '(a)') mvin
    read( filein, '(a)') mvout
endif

close( filein )

open(file=filename_box, status='old', newunit=filein)
read(filein,*) h(1,1), h(1,2), h(1,3)
read(filein,*) h(2,1), h(2,2), h(2,3)
read(filein,*) h(3,1), h(3,2), h(3,3)
close(filein)

open( file=filename, status='old', newunit=filein )
do i=1,nspec
    associate( this_spec => spec(i) )
        do j=1, this_spec%num
            associate( this_atom => this_spec%atom(j) ) 
                read(filein,*) rtemp
                this_atom%r = cell_to_lab( rtemp, h ) ! Convert from cell to lab. coordinates.
            end associate
        end do
        if (this_spec%read_flags ) then
            open( file=this_spec%flag_filename, status='old', newunit=this_spec%fileflags )
            read( this_spec%fileflags,* ) this_spec%atom%flag
        else
            this_spec%atom%flag = 1
        endif
    end associate
end do
close(filein)

open (file=filename_out,form='formatted', newunit=fileout)

! Write out 'atom....' details.
do j=1,nspec
    do i=1, spec(j)%num
        associate( this_atom => spec(j)%atom(i) )
        if (spec(j)%read_flags) then
            write(fileout,*)'atom    ', spec(j)%name, real( this_atom%r(1) ), real( this_atom%r(2) ), real( this_atom%r(3) ), this_atom%flag
        else
            write(fileout,*)'atom    ', spec(j)%name, real( this_atom%r(1) ), real( this_atom%r(2) ), real( this_atom%r(3) )
        end if
        end associate
    end do
end do

! Write out 'spec....' details.
do j=1, nspec
    write(fileout,901) spec(j)%name, real(spec(j)%radius), real(spec(j)%rgb(1)), real(spec(j)%rgb(2)), real(spec(j)%rgb(3))
end do

! Write out 'bond....' details.
do j=1, npairs
    write(fileout,*)'bonds    ', pair(j)%name1, pair(j)%name2, real(pair(j)%minl), real(pair(j)%maxl), real(pair(j)%bondradius), real(pair(j)%bondgrey)
end do
 
! Complete file.
write(fileout,*)'tmat -1.0 0.0 0.0 -0.0 -1.0 0.0 0.0 0.0 1.0' !
write(fileout,*)'dist    15.000'                              !
write(fileout,*)'inc      3.000'                              !
write(fileout,*)'scale   10.000'                              !
write(fileout,*)'rfac 1.00'                                   !
write(fileout,*)'bfac 2.00'                                   !
write(fileout,*)'gramp    0.000    0.000'                     !
write(fileout,*)'switches 1 0 1 0 0 1 1 0 0'                  !
    
992 format (3(f9.6,1x),'t',i4,'  c  ',A10,' r ',F9.6)
901 format ('spec', t10, a3, t17, f6.3, t25, f6.3, t32, f6.3, t42, f6.3)

close(fileout)

if (write_movie) then

    open( file=filename_box, status='old', newunit=filebox )
    open( file=mvin, status='old', newunit=filein )
    open( file=mvout, newunit=fileout )
    do i=1, nframes * nskip

        read( filebox, * ) h(1,1), h(1,2), h(1,3)
        read( filebox, * ) h(2,1), h(2,2), h(2,3)
        read( filebox, * ) h(3,1), h(3,2), h(3,3)
        read( filebox, * ) boxlen

        if (mod(i,nskip) == 0) then
            write(fileout,*)'frame        t', i
            do j=1, nspec
                associate( this_spec => spec(j) )
                    if (this_spec%read_flags ) read( this_spec%fileflags,* ) this_spec%atom%flag
                    do k=1, this_spec%num
                        associate( this_atom => spec(j)%atom(k) )
                            read(filein,*) rtemp
                            this_atom%r = cell_to_lab( rtemp, h )
                            if ( this_spec%read_flags ) then
                                write(fileout,*) real( this_atom%r(1) ), real( this_atom%r(2) ), real( this_atom%r(3) ), this_atom%flag
                            else
                                write(fileout,*) real( this_atom%r(1) ), real( this_atom%r(2) ), real( this_atom%r(3) )
                            end if
                        end associate
                    end do
                end associate
            end do
        else
            do j=1, nspec
                do k=1, spec(j)%num
                    read(filein,*) rtemp
                end do
            end do
        endif 
    end do
    close( filein )
    close( fileout )
    close( filebox )
endif

end program

function cell_to_lab( r, h )
    implicit none
    double precision, dimension(3), intent(in) :: r
    double precision, dimension(3,3), intent(in) :: h
    double precision, dimension(3) :: cell_to_lab
    integer :: i
    forall (i=1:3) cell_to_lab(i) = h(i,1)*r(1) + h(i,2)*r(2) + h(i,3)*r(3)
end function cell_to_lab
