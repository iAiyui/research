PROGRAM main
! disable implicit definition
implicit none

!定数項-モジュール化したい-
real*4,   parameter    :: dt =0.00004!0.000208
real*4,   parameter    :: dx = 0.02!0.1

real*4,parameter       :: tiv = 0.005!解析時間長
real*4,parameter    :: reg_lx =  6.0       !空間の横幅[m] real*4,parameter    :: reg_ly =  15.0       !空間の長さ[m]
real*4,parameter    :: reg_ly =  6.0       !空間の横幅[m] real*4,parameter    :: reg_ly =  15.0       !空間の長さ[m]
real*4,parameter       :: frequency = 25.4       !純音の場合の印加する周波数

real*4,parameter    :: reg_ix = reg_lx/dx       !教室セル数[n]
real*4,parameter    :: reg_jx = reg_ly/dx       !教室セル数[n]
real*4,parameter    :: span = 0.5                   !受音点ピッチ[m]

integer*4,parameter	:: ix = (reg_lx)/dx         !401! 1000number of x-directional spatial intervals
integer*4,parameter	:: jx = (reg_ly)/dx     !901number of y-directional spatial intervals
integer*4,parameter	:: tx = (tiv/dt)                   !計算秒数[s] 以前は更新回数[n]

real*4,parameter	:: td = (1/frequency)/dt			! input frequency 純音の時使います
integer*4,parameter	:: id = (ix)/2    !左右中心 x-directional position of input
integer*4,parameter	:: jd = (jx)/2!で設定音源		   !奥の壁から7m y-directional position of input

real*4,	 parameter	:: crn = 0.7		! Courant number (= phase velocity * time interval / spatial interval)
real*4,   parameter	:: dd = 199.0/200.0	! parameter of Higdon's absorption boundary


real*4, 	dimension(0:ix+1,0:jx+1)	:: p1	             ! sound pressure (t = n + 1 / 2)
real*4, 	dimension(0:ix+1,0:jx+1)	:: p2	             ! sound pressure (t = n - 1 / 2)
real*4, 	dimension(0:ix+1,0:jx+1)	:: p3	             ! sound pressure (t = n - 3 / 2)
real*4, 	dimension(0:ix+1,0:jx+1)	:: u1	             ! x-directional velocity (t = n + 1)
real*4, 	dimension(0:ix+1,0:jx+1)	:: u2	             ! x-directional velocity (t = n)
real*4, 	dimension(0:ix+1,0:jx+1)	:: v1	             ! y-directional velocity (t = n + 1)
real*4, 	dimension(0:ix+1,0:jx+1)	:: v2	             ! y-directional velocity (t = n)
real*4, 	dimension(tx)				:: pin               ! input wave
real*4,   dimension(2)         ::updw
real*4,   dimension(2)         :: xlim
!動的配列の確保
!座標情報(x座標，y座標)
integer*4, allocatable, dimension(:)   :: x                   ! coordinate_x
integer*4, allocatable, dimension(:)   :: y                   ! coordinate_y
!
!!直線のベクトル情報(各軸方向ベクトル，単位ベクトル，判定用変数)
!real*4,    allocatable, dimension(:)   ::vectorx              ! vector between 2coordinate
!real*4,    allocatable, dimension(:)   ::vectory
!real*4,    allocatable, dimension(:)   ::vectorSx
!real*4,    allocatable, dimension(:)   ::vectorSy
!real*4,    allocatable, dimension(:)   ::Judge      
!
!!直線情報(傾き，高さ,長さ)
real*4,    allocatable, dimension(:)   ::katamuki_x           !slope between the 2points
real*4,    allocatable, dimension(:)   ::katamuki_y           !slope between the 2points
real*4,    allocatable, dimension(:)   ::katamuki             !slope between the 2points
real*4,    allocatable, dimension(:)   ::high                 !slope between the 2points
real*4,    allocatable, dimension(:)   :: norm                ! vector_Absolute
!
!!挟み込む範囲
!integer*4, allocatable, dimension(:)   :: x2                  ! coordinate_x
!integer*4, allocatable, dimension(:)   :: y2 
real*4,    allocatable, dimension(:,:) ::tmp
real*4,    allocatable, dimension(:)   ::tmpno
!
!!直線座標算出用
real*4,    allocatable, dimension(:,:,:) ::length
!real*4,    allocatable, dimension(:,:)   :: cor_x             ! insert the number of y axis cells[m]
!real*4,    allocatable, dimension(:,:)   :: cor_y             ! insert the number of x axis cells[m]
!                                    
!integer*4, allocatable, dimension(:,:)   :: cor
!integer*4, allocatable, dimension(:,:)   :: cor_x_round       ! insert the number of y axis cells
!integer*4, allocatable, dimension(:,:)   :: cor_y_round       ! insert the number of x axis cells
!
!integer*4, allocatable, dimension(:,:)   :: range_x_max
!integer*4, allocatable, dimension(:,:)   :: range_x_min
!
!
real*4,    allocatable, dimension(:)   :: unit_vector_x       !unit vector between the 2points
real*4,    allocatable, dimension(:)   :: unit_vector_y
real*4,    allocatable, dimension(:)   :: nx
real*4,    allocatable, dimension(:)   :: ny   
!
!
!
!
character*256	:: fn, str                  ! output file name and temporal pass string
!character*20    :: filename					!ファイル名書き出し
!integer*4  :: sampleno
integer*4 t1, t2, t_rate, t_max, diff      ! time measurement

integer*4   :: coordinateno, no 
integer*4	:: t, i, j, l, k, counter 
integer*4   :: pos_y, pos_x, filenumber  	! loop variable
integer*4 :: source , output                        !which source to use 
real*4		:: pai							             ! circular constant
real*4		:: a, b, c, d, e			         	! coefficients of Higdon's absorption boundary
real*4    :: Z , arufa , a1                !cal abs
real*4    :: xx,yy,r,Rlength               !cal src
!------------------------------------------------------------+
!reading coodinate
!------------------------------------------------------------+
open (50,file='C:\Users\n\Documents\16_program\input\4_cor.txt')
    read (50,'(I3)')  coordinateno
    allocate( x       (coordinateno+1))
    allocate( y       (coordinateno+1))
    
!    allocate( x2       (coordinateno))
!    allocate( y2       (coordinateno))
    
    allocate( length        (coordinateno,ix+1,jx+1) )
    allocate( tmp           (coordinateno,ix+1) )
    allocate( tmpno         (coordinateno) )
!    allocate( cor_x         (coordinateno,ix+1))      !座標数=直線の本数の関係を使用. (x,y) = (直線の本数,x方向への空間広さ)
!    allocate( cor_y         (coordinateno,jx+1))       !ix,jx配列長にふさわしい方を判定する仕組みが必要
!    allocate( cor_x_round   (coordinateno,ix+1))
!    allocate( cor_y_round   (coordinateno,jx+1))
!    allocate( cor           (coordinateno,ix+1))

    allocate( norm            (coordinateno))
    allocate( high            (coordinateno))

 !   allocate( vectorx         (coordinateno))
 !   allocate( vectory         (coordinateno))
 !   allocate( vectorSx        (coordinateno))
 !   allocate( vectorSy        (coordinateno))
 !   allocate( Judge           (coordinateno))

    allocate( unit_vector_x   (coordinateno))
    allocate( unit_vector_y   (coordinateno))
    allocate( nx              (coordinateno))    
    allocate( ny              (coordinateno))    

    allocate( katamuki_x      (coordinateno))
    allocate( katamuki_y      (coordinateno))
    allocate( katamuki        (coordinateno))
    

 !   allocate( range_x_max     (0:ix+1))
 !   allocate( range_x_min     (0:ix+1))
    
    rewind(50)
    print*, coordinateno
    !--------------------------------------------------------------------------------+
    !以下の'10I3'となっている部分を10->13とか読み込める範囲広げれば多角形領域作れます/0215yui
    !--------------------------------------------------------------------------------+
    read(50,'(3x,10I3)') (x(i), y(i), i=1,coordinateno+1) 
    !print*,x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),x(5),y(5),coordinateno
close(50)



!integer*4 ::grid                           !cal rp axis

!data rp_x / 25/    !受音点ピッチ(x方向)
!data rp_y / 25/    !受音点ピッチ(y方向)
data counter/20/   !生成ファイル数(装置番号)
data arufa/0.8/    !壁面吸音率1まで
data source/1 /    !0:純音 1:インパルス
data Rlength/7.65/ !インパルス音,周波数特性の決定 7.65:500Hz, 1.91:2k, 9.61:0.1グリッド@500Hz
data output / 1 /  !0=波形のみ　1 =伝搬出力

call cal_abs(arufa, Z) 


! indicate output file name
fn	= 'sound-wave-0'

! constants
pai	= acos(-1.0)
a	= 2.0 * (crn - 1.0) / (crn + 1.0)
b	= ((crn - 1.0) / (crn + 1.0)) ** 2
c	= 2.0 * (crn - 1.0) / (crn + 1.0) * dd
d	= dd * 2.0
e	= dd ** 2

! initialize
do i = 0, ix + 1
do j = 0, jx + 1
	p1(i,j)	= 0.0
	p2(i,j)	= 0.0
	p3(i,j)	= 0.0
	u1(i,j)	= 0.0
	u2(i,j)	= 0.0
	v1(i,j)	= 0.0
	v2(i,j)	= 0.0
end do
end do

!do i = 1,coordinateno
!    do j = 1,ix+1
!        cor_x(i,j) = 0
!        cor_x_round(i,j) =1
!    end do
!
!    do j = 1,jx+1
!        cor_y(i,j) = 0
!        cor_y_round(i,j) =1
!    end do
!    do j = 1,ix
!        cor(i,j) =0
!    end do
!end do
!
!do i = 1,coordinateno
!    do j = 1,ix+1
!        cor_x(i,j) = 0
!        cor_x_round(i,j) =1
!    end do
!
!    do j = 1,jx+1
!        cor_y(i,j) = 0
!        cor_y_round(i,j) =1
!    end do
!    do j = 1,ix
!        cor(i,j) =0
!    end do
!end do

do i = 1,coordinateno
    do j= 1,ix+1
        do k = 1,jx+1
        length(i,j,k) = 0
        end do
    end do
end do

do k = 1,coordinateno
    do i = 1, ix+1
        tmp(k,i) = 0
    end do
end do


!------------------------------------------------------------+
! input wave definition (one cycle of sinusoidal wave on one point)
!------------------------------------------------------------+
!First source
if(source == 0)then
    do t = 1,tx! 2*td!一波長だけ印可するtx
    	if(t <= tx)then!2*td) then
    		pin(t)	= sin(real(t) / real(td) * 2.0 * pai) 
      else
        pin(t)	= 0.0
    	end if
    end do
else if(source == 1)then
    !Second sourc
    do i = 2,ix
    do j = 2,jx
        xx = (i-id)
        yy = (j-jd)
        r = sqrt(xx**2+yy**2)
        p2(i,j) = exp(-(r/Rlength)**2)
        end do
    end do
end if


!前準備
!-----------------------------------------+
!get slope between 2points
!-----------------------------------------+
!t =1
do t = 1,coordinateno
         if (x(t) /= x(t+1) .and. y(t) /= y(t+1) )then 
             call cal_katamuki(x(t), x(t+1), y(t), y(t+1), high(t), norm(t), katamuki_x(t), katamuki_y(t), katamuki(t) )
         else if(x(t) == x(t+1) )then
              call cal_x_same(x(t), x(t+1), y(t), y(t+1), high(t), norm(t), katamuki_x(t), katamuki_y(t), katamuki(t) )
         else if(y(t) == y(t+1) )then
              call cal_y_same(x(t), x(t+1), y(t), y(t+1), high(t), norm(t), katamuki_x(t), katamuki_y(t), katamuki(t) )
         end if
    !get Normal vector by turn -90degree unit vector 単位ベクトル(t)から-90回した法線ベクトル取得
    call cal_unit_vector  ( katamuki_x(t), katamuki_y(t), norm(t), unit_vector_x(t), unit_vector_y(t) )
    call cal_linertransform(unit_vector_x(t), unit_vector_y(t), nx(t), ny(t) )
end do

!-----------------------------------------+
!remake:直線のサンプリング!
!       配列要素を1とそれ以外に分ける
!-----------------------------------------+

do i = 1,coordinateno 
    if (abs(katamuki(i) ) >1 .and. x(i) /= x(i+1) )then
        if( y(i) < y(i+1) )then
            do pos_y = y(i),y(i+1)
                length(i, int( (pos_y - high(i))/katamuki(i)), pos_y ) = 1
            end do
        else
            do pos_y = y(i), y(i+1), -1
                length(i, int( (pos_y - high(i))/katamuki(i)), pos_y ) = 1
            end do
        end if
    
    else if( abs(katamuki(i) ) < 1)then
        if( x(i) < x(i+1) )then
            do pos_x = x(i), x(i+1)
                length(i, pos_x, int( katamuki(i)*pos_x + high(i) )) = 1
            end do
        else
            do pos_x = x(i),x(i+1),-1
                length(i, pos_x, int(katamuki(i)*pos_x + high(i) )) = 1
            end do
        end if
    
    else if(x(i) == x(i+1) )then
        if( y(i) < y(i+1) )then
            do pos_y = y(i),y(i+1)
                length(i, x(i), pos_y ) = 1
            end do
        else
            do pos_y = y(i), y(i+1), -1
                length(i, x(i), pos_y ) = 1
            end do
        end if
    end if
end do




do k =1,coordinateno
do i = 1,ix+1
do j = 1,jx+1
if (length(k, i, j ) == 1)then
    call gen_rp(20, 1,'a',1,real(k),real(i),real(j) )
    !print *, k, i, j
end if
end do
end do
end do



!-----------------------------------------+
!  old  ソートして挟み込む範囲の決定■■■!
!-----------------------------------------+
    do i = 1,ix
    counter = 1
        do k =1,coordinateno
            do j = 1,jx
                if(length(k,i,j) == 1)then
                    tmp(counter ,i) = j
                    print *, counter,i, tmp(counter, i) 
                    counter = counter +1
                end if
            end do
        end do
    end do
print *, tmp(1,161),tmp(2,161),tmp(3,161),tmp(4,161)
print *, tmp(1,155),tmp(2,155),tmp(3,155),tmp(4,155)

    do i = 1,ix-1
        do no = 1, coordinateno
            tmpno(no) = tmp(no, i)
        end do
        call countup_no(tmpno, coordinateno, updw,2)
        print*, updw(1),',', updw(2)
    end do


!-----------------------------------------+
!■■■x座標をソートして挟み込む範囲の決定■■■!
!-----------------------------------------+
call coo_maxmin(x, coordinateno, xlim)


call system_clock(t1)   ! 開始時を記録
!----------------------------------------+
! time loop
!----------------------------------------+
do t = 1, tx
! update of sound pressure
    do i = xlim(1),xlim(2)
        do j = 
            p1(i,j)	= p2(i,j) - crn * (u2(i,j) - u2(i-1,j) + v2(i,j) - v2(i,j-1))
        end do
    end do
    
    !input of sound pressure
    !インパルス応答使うときはコメントアウト
    if(source == 0)then
        p1(id,jd) = pin(t) + p1(id,jd)
    end if
    
    ! swap of sound pressure
    do i = xlim(1)-1,xlim(2)+1!0, ix + 1
        do j =
            p3(i,j)	= p2(i,j)
            p2(i,j)	= p1(i,j)
        end do
    end do
    
    
    ! update of x-directional velocity
    do i = xlim(1)-1,xlim(2)+1
        do j =
            u1(i,j)	= u2(i,j) - crn * (p2(i+1,j) - p2(i,j))
        end do
    end do
    
    ! update of y-directional velocity
    do i = xlim(1),xlim(2)!1, ix
        do j =
            v1(i,j)	= v2(i,j) - crn * (p2(i,j+1) - p2(i,j))
        end do
    end do





end program main





!■■■■■■■■吸音率の計算■■■■■■■■■■■■■■■■■■■■■!
!in:吸音率arufa
!
!out:インピーダンス（エネルギー量)Z
!
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine cal_abs( arufa, Z)
implicit none 

real*4,intent(in)  :: arufa
real*4,intent(out) :: Z

real tmp
tmp=1-arufa
Z = ((1+SQRT(tmp))/(1-SQRT(tmp)))

return 
end subroutine 


!■■■■■■■■2点間の傾き,絶対値を求める■■■■■■■■■■■■■■■■■■!
!in
!2組の座標 A( x(1),y(1) ),B( x(2),y(2) )
!out
!・2点間の距離(絶対値) norm
!・ベクトルxの長さ,ベクトルyの長さ katamuki_x, katamuki_y
!・傾き katamuki
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine cal_katamuki (x1, x2, y1, y2, high, norm, katamuki_x, katamuki_y, katamuki )
implicit none

integer*4 x1, x2, y1, y2
real*4 norm
real*4 katamuki_x, katamuki_y, katamuki
real*4 high

katamuki_x = x2 -x1
katamuki_y = y2 -y1
norm = sqrt( (katamuki_x )**2 +( katamuki_y )**2 )
katamuki   = real(katamuki_y/katamuki_x)

x1 = int(x1)
y1 = int(y1)
call cal_high(x1, y1, katamuki, high) !高さを出力
return 
end subroutine

!■■■■■■■x(t) = x(t+1)のときに値を代入する■■■■■■■■■■■■■■■■!
!in
!2組の座標 A( x(1),y(1) ),B( x(2),y(2) )
!out
!・2点間の距離(絶対値) norm
!・ベクトルxの長さ,ベクトルyの長さ katamuki_x, katamuki_y
!・傾き katamuki
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine cal_y_same (x1, x2, y1, y2, high, norm, katamuki_x, katamuki_y, katamuki )
implicit none

integer*4 x1, x2, y1, y2
real*4 norm
real*4 katamuki_x, katamuki_y, katamuki
real*4 high


katamuki = 0
katamuki_y =1
katamuki_x = abs( x2 -x1 )
norm = katamuki_x
    
high = y1

x1 = int(x1)
y1 = int(y1)
return 
end subroutine


!■■■■■■■x(t) = x(t+1)のときに値を代入する■■■■■■■■■■■■■■■■!
!in
!2組の座標 A( x(1),y(1) ),B( x(2),y(2) )
!out
!・2点間の距離(絶対値) norm
!・ベクトルxの長さ,ベクトルyの長さ katamuki_x, katamuki_y
!・傾き katamuki
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine cal_x_same (x1, x2, y1, y2, high, norm, katamuki_x, katamuki_y, katamuki )
implicit none

integer*4 x1, x2, y1, y2
real*4 norm
real*4 katamuki_x, katamuki_y, katamuki
real*4 high


katamuki = 1000
katamuki_x = 0
katamuki_y = abs( y2 - y1 )
norm = katamuki_y

high = 1000000


x1 = int(x1)
y1 = int(y1)
return 
end subroutine

!■■■■■■■■2点を結ぶ直線の高さを求める■■■■■■■■■■■■■■■■■■!
!in
!1組の座標 A( x(1),y(1) ),傾きkatamuki
!out
!直線の高さhigh
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine cal_high(x1, y1, slp, high)

integer*4,intent(in) ::x1, y1
real*4,intent(in) ::slp
real*4,intent(out):: high

high = real(-x1*slp+y1)

return
end subroutine

!■■■■■■■■単位ベクトルを求める■■■■■■■■■■■■■■■■■■■■■!
!in
!katamuki_x,katamuki_y , normAB
!out             
!・単位ベクトルunit_vector_x, unit_vector_y
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine cal_unit_vector(katamuki_x, katamuki_y, norm, unit_vector_x, unit_vector_y)

real*4 unit_vector_x, unit_vector_y
real*4 katamuki_x, katamuki_y
real*4 norm

unit_vector_x = katamuki_x/norm
unit_vector_y = katamuki_y/norm

return
end

!■■■単位ベクトルを-90度方向へ回転させた法線ベクトルを求める■■■■■■■■■!
!in
!unit_vector_x,unit_vector_y
!out
!・法線ベクトルnx,ny
!
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine cal_linertransform(unit_vector_x,unit_vector_y,nx,ny)

real*4 unit_vector_x,unit_vector_y,nx,ny

nx = ( 0 * unit_vector_x  + 1 * unit_vector_y) !ここあってる？
ny = ( -1 * unit_vector_x + 0 * unit_vector_y)

return
end 

!■■■降順ソート■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
!in
!配列長N, 配列要素data
!out
!並び替えた配列要素dataを割り当てて返す
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine sortdp(N,data)
  implicit none
  integer::i,j,N
  integer*4::tmp,data(1:N)

  do i=1,N-1
     do j=i+1,N
        if(data(i) > data(j))then
           tmp=data(i)
           data(i)=data(j)
           data(j)=tmp
        end if
     end do
  end do

  return
end subroutine 

!■■■昇順ソート■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
!in
!配列長N, 配列要素data
!out
!並び替えた配列要素dataを割り当てて返す
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine sortao(N,data)
  implicit none
  integer::i,j,N
  integer*4::tmp,data(1:N)

  do i=1,N-1
     do j=i+1,N
        if(data(i) < data(j))then
           tmp=data(i)
           data(i)=data(j)
           data(j)=tmp
        end if
     end do
  end do

  return
end subroutine 

!■■■降順ソート■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
!in
!配列長N, 配列要素data
!out
!配列中の0以外の要素数ele
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine count_N(N, data, ele)
    implicit none

    integer :: i, N
    integer*4, intent(in) :: data(1:N)
    integer*4,intent(out):: ele

    ele =0
    do i = 1,4
        if(data(i) /= 0)then
            ele = ele +1
        end if
    end do
    return
end subroutine count_N

!■■書き出し■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
!in
!rp,sts, wrt
!out
!rp:受音点数
!sts:どういう風に生成したファイルを返すのか(1:開きっぱなし,2:閉じて返す)←必要？
!wrt:書き込む内容
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine gen_rp(rp, sts, filename,lng,  wrt,wrt2,wrt3)
    implicit none
    integer*4,intent(in) :: rp,sts, lng
    real*4   ,intent(in) :: wrt,wrt2,wrt3
    character(lng),intent(in):: filename

    open(rp, file = 'data_imp\' //trim(filename)// '',position = 'APPEND')
    write(rp,*) wrt,',',wrt2,',',wrt3
    close(rp)

    return
end subroutine 
!-----------------------------------------------------------------------------+
!配列の最大値と最小値のみ抽出して返すサブルーチン
!in:list
!
!out:list[1]へ最大値，llist[2]へ最小値，あとは0で埋めた元の配列
!------------------------------------------------------------------------------+
subroutine list_min_max(list, N)
    implicit none
    integer*4 :: N,i
    integer*4,intent(inout) ::list(1:N)

    list(1) = maxval(list)
    list(2) = minval(list)
    
    do i =3,N
        list(i) = 0
    end do

    return
end subroutine list_min_max


!-------------------------------------------------------+
!
!-------------------------------------------------------+
subroutine countup_no(tmp1, N, updw, no)
    implicit none
    integer*4              ::N,i,counter,no
    real*4,intent(in)      ::tmp1(1:N)
    real*4,intent(inout)   ::updw(no)


    counter = 0
    do i =1,N
        if(tmp1(i) /= 0)then
            counter = counter +1
        end if
    end do
    call lid_zero(updw, tmp1, counter,N )
    return
end subroutine

subroutine lid_zero(updw, list, lng, N )
    implicit none
    integer*4         :: lng, N, i, counter
    real*4,intent(in) :: list(1:N)  !元の配列0含んでる
    integer*4         :: tmp2(1:lng)        !この関数内でのみ使う配列
    real*4,intent(inout)    :: updw(2)
    
    counter = 0
    do i =1,N
        if(list(i) /= 0)then
            counter =counter +1 
            tmp2(counter) = list(i)
        end if
    end do
    call coo_maxmin(tmp2, lng, updw)
    return
end subroutine


subroutine coo_maxmin(list, lng, updw)
    implicit none
    integer*4    :: lng,i
    integer*4    :: list(1:lng)
    real*4,intent(inout)    :: updw(2)
    
    updw(1) = minval(list)
    updw(2) = maxval(list)
    return
end subroutine

!■■■書き出し内容整形■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
!in
!rp,sts, wrt
!out
!rp:受音点数
!sts:どういう風に生成したファイルを返すのか(1:開きっぱなし,2:閉じて返す)
!wrt:書き込む内容
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
!subroutine 
!    implicit none
!    integer*4,intent(in) :: rp,sts, lng
!    real*4   ,intent(in) :: wrt
!    character(lng),intent(in):: filename
!
!    open(rp, file = 'data_imp\' //trim(filename)// '.txt',position = 'APPEND')
!    write(rp,'(f6.2)') wrt
!    close(rp)
!
!    return
!end subroutine gen_rp
