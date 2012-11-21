!===============================================================================
!
!  plate_motion is a set of F90 functions and subroutines for calculating relative
!  and absolute plate motion on a sphere.
!
!  Included are NUVEL-1 and NUVEL-1A, as well GSRM plate velocities.
!
!  Andy Nowacki, University of Bristol, 2011/09
!
!-------------------------------------------------------------------------------
!  History:
!     + 2011/11/09:  Added NNR-MORVEL(±56) plate velocities.
!     + 2011/11/09:  Beginning of support for plate boundaries.  This will allow
!                    for the auto-determination of which plate a given point is in.
!===============================================================================
module plate_motion
!===============================================================================

      implicit none
      
!  ** Hide the internal utility toupper function
      private :: toupper

!  ** size constants
      integer, parameter, private :: i4 = selected_int_kind(9) ; ! long int
      integer, parameter, private :: r4 = selected_real_kind(6,37) ; ! SP
      integer, parameter, private :: r8 = selected_real_kind(15,307) ; ! DP
      
!  ** precision selector
      integer, parameter, private :: rs = r8
      
!  ** maths constants and other useful things
      real(rs), parameter, private :: pi = 3.141592653589793238462643 ;
      real(rs), parameter, private :: to_rad = 1.74532925199433e-002 ;  
      real(rs), parameter, private :: to_deg = 57.2957795130823e0 ;  
      real(rs), parameter, private :: to_km = 111.194926644559 ;      
      real(rs), parameter, private :: big_number = 10.e36 ;     
      real(rs), parameter, private :: R_EARTH = 6371._rs
      
!  ** IO
      integer, parameter, private :: lu_stdin = 5
      integer, parameter, private :: lu_stdout = 6
      integer, parameter, private :: lu_stderr = 0

!  ** Parameters for the module
      integer, parameter, public :: MODEL_NAME_LENGTH = 12
      integer, parameter, public :: PLATE_NAME_LENGTH = 2
      integer, parameter, public :: LONG_PLATE_NAME_LENGTH = 30
      integer, parameter, private :: MAXIMUM_PLATES = 100
      integer, parameter, private :: MAXIMUM_VERTICES = 2000

!  ** Available models in the plate_motion module
      integer, parameter, private :: N_AVAILABLE_MODELS = 5
      integer, parameter, private :: N_PLATES = 56
            
!  ** Derived types
!  Plate type, with abbreviation, Euler pole and rotation rate
      type plate
         character(len=PLATE_NAME_LENGTH) :: name  ! Short plate name
!         character(len=LONG_PLATE_NAME_LENGTH) :: longname ! Long name
         real(rs) :: pole(3)      ! Lon, lat, rotation rate (deg/Ma)of Euler pole 
!         real(rs) :: vertices(MAXIMUM_VERTICES,2) ! lon,lat of plate boundaries
      end type plate
      
!  APM model type, containing several plates, its name and others
      type apm_model
         character(len=MODEL_NAME_LENGTH) :: name    ! Model name
         integer :: n_plates                         ! Number of plates in model
         type(plate), dimension(MAXIMUM_PLATES) :: p ! Plate of type(plate)
      end type
      
!  List of plate names, taken from the NNR-MORVEL56 list
      type plate_names_list
         character(len=PLATE_NAME_LENGTH) :: ab(N_PLATES)
         character(len=LONG_PLATE_NAME_LENGTH) :: name(N_PLATES)
      end type
      character(len=LONG_PLATE_NAME_LENGTH), parameter :: NO_PLATE_NAME = '--'
      
      
!===============================================================================
!  PLATE DEFINITIONS
!  The models are declared as a parameter visible to all functions and subroutines
!  in the module; so each is an element of an array of the apm_model derived type.
!  Declare the models here using the type declaration syntax: 
!        x = type1(type2(element1, element2), (/ array(1), array(2) /), 'character')
!  Follow the brackets carefully with this syntax and don't forget the &s!
!  Units for poles are lon (deg), lat (deg), rate (deg/Ma).
!
!  REFERENCES:
!     HS2-NUVEL1:   Gripp & Gordon, 1990. Current plate velocities relative to
!                     the hostpots incorporating the NUVEL-1 global plate motion
!                     model. GRL, 17, 1109-1112.
!     HS3-NUVEL1A:  Gripp & Gordon, 2002. Young tracks of hostpots and current
!                     plate velocities. GJI, 150, 321-361.
!     GSRM-NNRv1.0: Creemer, Holt and Haines, 2003, An integrated global model 
!                     of present-day plate motions and plate boundary deformation.
!                     GJI, 154, 8-34.
!     NNR-MORVEL56: Argus, Gordon & DeMets, 2011. Geologically current motion
!     & NNR-MORVEL    of 56 plates relative to the no-net-rotation reference
!                     frame.  G^3, 12, Q11001.
!===============================================================================
      type(apm_model), parameter, private, dimension(N_AVAILABLE_MODELS) :: &
MODEL = (/                                              &
   apm_model (                                          &
      'HS2-NUVEL1',                                     &
      14,                                               &
      (/                                                &
         plate( "AF", (/    3.600,  -5.500, 0.1500 /)), &
         plate( "AN", (/   65.900, -14.800, 0.1100 /)), &
         plate( "AR", (/   18.400,  16.800, 0.5400 /)), &
         plate( "AU", (/   41.800,   9.600, 0.7600 /)), &
         plate( "CA", (/   -5.700, -62.400, 0.1700 /)), &
         plate( "CO", (/ -115.900,  18.400, 1.2900 /)), &
         plate( "EU", (/   58.100, -44.800, 0.0900 /)), &
         plate( "IN", (/   21.900,  16.600, 0.5500 /)), &
         plate( "JF", (/   60.000, -34.800, 0.9500 /)), &
         plate( "NA", (/  -11.100, -67.200, 0.2800 /)), &
         plate( "NZ", (/  -90.200,  45.700, 0.4600 /)), &
         plate( "PA", (/   90.000, -60.200, 0.9800 /)), &
         plate( "PH", (/  -19.900, -49.400, 1.1100 /)), &
         plate( "SA", (/   74.700, -70.300, 0.3200 /))  &
         /)                                             &
   ),                                                   &
&                                                       &
   apm_model (                                          &
      'HS3-NUVEL1A',                                    &
      15,                                               &
      (/                                                &
         plate( "AF", (/   21.136, -43.386, 0.1987 /)), &
         plate( "AN", (/   74.514, -47.339, 0.2024 /)), &
         plate( "AR", (/   23.175,   2.951, 0.5083 /)), &
         plate( "AU", (/   44.482,  -0.091, 0.7467 /)), &
         plate( "CA", (/   25.925, -73.212, 0.2827 /)), &
         plate( "CO", (/ -116.997,  13.171, 1.1621 /)), &
         plate( "EU", (/   73.474, -61.901, 0.2047 /)), &
         plate( "IN", (/   26.467,   3.069, 0.5211 /)), &
         plate( "JF", (/   61.633, -39.211, 1.0122 /)), &
         plate( "NZ", (/  -90.913,  35.879, 0.3231 /)), &
         plate( "NA", (/   13.400, -74.705, 0.3835 /)), &
         plate( "PA", (/   90.326, -61.467, 1.0613 /)), &
         plate( "PH", (/  -16.668, -53.880, 1.1543 /)), &
         plate( "SS", (/   52.228, -76.912, 0.4451 /)), &
         plate( "SA", (/   80.401, -70.583, 0.4358 /))  &
      /)                                                &
   ),                                                   &
&                                                       &
   apm_model (                                          &
      'GSRM-NNRv1.0',                                   &
      24,                                               &
      (/                                                &
         plate( "AM", (/ -103.800,  60.000, 0.3020 /)), &
         plate( "AN", (/ -120.600,  62.700, 0.2420 /)), &
         plate( "AR", (/   -2.700,  50.600, 0.5500 /)), &
         plate( "AT", (/   27.000,  42.000, 1.4110 /)), &
         plate( "AU", (/   35.800,  34.200, 0.6260 /)), &
         plate( "CA", (/  -91.100,  37.000, 0.3020 /)), &
         plate( "CO", (/  -63.683,  25.558, 1.5024 /)), &
         plate( "CR", (/  -15.790, -56.338, 0.8964 /)), &
         plate( "EU", (/  -97.400,  56.400, 0.2790 /)), &
         plate( "IN", (/  -13.700,  52.600, 0.4890 /)), &
         plate( "JF", (/   54.043, -36.904, 0.8843 /)), &
         plate( "NA", (/  -82.300,   1.700, 0.2110 /)), &
         plate( "NU", (/  -78.300,  52.600, 0.2940 /)), &
         plate( "NZ", (/  -98.400,  44.400, 0.6510 /)), &
         plate( "OK", (/  -43.900, -38.700, 1.0370 /)), &
         plate( "PA", (/  105.500, -64.900, 0.6500 /)), &
         plate( "PH", (/  -28.900, -45.100, 0.8740 /)), &
         plate( "RV", (/  -73.089,  20.546, 4.4808 /)), &
         plate( "SA", (/ -119.500, -14.500, 0.1140 /)), &
         plate( "SC", (/ -163.100,  73.800, 0.4110 /)), &
         plate( "SM", (/  -63.300,  57.300, 0.2820 /)), &
         plate( "SN", (/  -90.200,  47.300, 0.3920 /)), &
         plate( "SS", (/  -77.676, -28.037, 0.1692 /)), &
         plate( "TA", (/  -86.800, -12.600, 0.6410 /))  &
      /)                                                &
   ),                                                   &
&                                                         &
   apm_model (                                            &
      'NNR-MORVEL',                                       &
      25,                                                 &
      (/                                                  &
         plate( "pa", (/  114.5565, -63.5635, 0.6497 /)), &
         plate( "am", (/ -122.5179,  63.3321, 0.2981 /)), &
         plate( "an", (/ -117.6938,  65.5957, 0.2509 /)), &
         plate( "ar", (/   -8.4816,  48.8922, 0.5602 /)), &
         plate( "au", (/   37.8817,  33.9135, 0.6327 /)), &
         plate( "ca", (/  -92.4244,  35.3723, 0.2869 /)), &
         plate( "co", (/ -124.2695,  26.9904, 1.1979 /)), &
         plate( "cp", (/   23.0387,  44.4615, 0.6093 /)), &
         plate( "eu", (/ -106.1893,  49.0695, 0.2234 /)), &
         plate( "in", (/   -3.2922,  50.3801, 0.5452 /)), &
         plate( "jf", (/   59.9816, -38.2434, 0.9509 /)), &
         plate( "lw", (/  -69.2798,  51.9679, 0.2867 /)), &
         plate( "na", (/  -80.4270,  -4.5494, 0.2088 /)), &
         plate( "nb", (/  -68.2247,  47.7696, 0.2931 /)), &
         plate( "mq", (/   11.0359,  49.1953, 1.1453 /)), &
         plate( "nz", (/ -100.9601,  46.3031, 0.6964 /)), &
         plate( "ps", (/  -31.3279, -45.9432, 0.9095 /)), &
         plate( "ri", (/ -107.2754,  20.2588, 4.5361 /)), &
         plate( "sa", (/ -112.3922, -22.1328, 0.1084 /)), &
         plate( "sc", (/ -105.8106,  22.9461, 0.1467 /)), &
         plate( "sm", (/  -84.3051,  50.0518, 0.3402 /)), &
         plate( "sr", (/ -110.8282, -32.0641, 0.1065 /)), &
         plate( "su", (/  -94.8062,  50.1768, 0.3377 /)), &
         plate( "sw", (/  -36.8459, -29.8875, 1.3616 /)), &
         plate( "yz", (/ -116.3321,  63.1609, 0.3343 /))  &
      /)                                                  &
   ),                                                     &
&                                                          &
   apm_model (                                             &
      'NNR-MORVEL56',                                      &
      56,                                                  &
      (/                                                   &
         plate( "pa", (/  114.6975, -63.5756,  0.6509 /)), &
         plate( "am", (/ -122.8242,  63.1704,  0.2973 /)), &
         plate( "an", (/ -118.1053,  65.4235,  0.2500 /)), &
         plate( "ar", (/   -8.4909,  48.8807,  0.5588 /)), &
         plate( "au", (/   37.9414,  33.8612,  0.6316 /)), &
         plate( "ca", (/  -92.6236,  35.1956,  0.2862 /)), &
         plate( "co", (/ -124.3074,  26.9346,  1.1978 /)), &
         plate( "cp", (/   23.0880,  44.4352,  0.6080 /)), &
         plate( "eu", (/ -106.5007,  48.8509,  0.2227 /)), &
         plate( "in", (/   -3.2898,  50.3722,  0.5438 /)), &
         plate( "jf", (/   60.0379, -38.3086,  0.9513 /)), &
         plate( "lw", (/  -69.5195,  51.8860,  0.2856 /)), &
         plate( "na", (/  -80.6447,  -4.8548,  0.2087 /)), &
         plate( "nb", (/  -68.4377,  47.6763,  0.2921 /)), &
         plate( "mq", (/   11.0524,  49.1891,  1.1440 /)), &
         plate( "nz", (/ -101.0564,  46.2348,  0.6957 /)), &
         plate( "ps", (/  -31.3615, -46.0242,  0.9098 /)), &
         plate( "ri", (/ -107.2861,  20.2450,  4.5359 /)), &
         plate( "sa", (/ -112.8327, -22.6179,  0.1090 /)), &
         plate( "sc", (/ -106.1485,  22.5244,  0.1464 /)), &
         plate( "sm", (/  -84.5154,  49.9506,  0.3393 /)), &
         plate( "sr", (/ -111.3224, -32.4957,  0.1072 /)), &
         plate( "su", (/  -95.0218,  50.0558,  0.3368 /)), &
         plate( "sw", (/  -36.8671, -29.9420,  1.3616 /)), &
         plate( "yz", (/ -116.6180,  63.0285,  0.3335 /)), &
         plate( "SL", (/ -143.4675,  50.7058,  0.2677 /)), &
         plate( "BH", (/  100.4994, -39.9983,  0.7988 /)), &
         plate( "MO", (/   92.6656,  14.2480,  0.7742 /)), &
         plate( "SS", (/  130.6236,  -2.8685,  1.7029 /)), &
         plate( "WL", (/  128.5186,   0.1050,  1.7444 /)), &
         plate( "CR", (/  170.5303, -20.3985,  3.9232 /)), &
         plate( "FT", (/  178.0679, -16.3322,  5.1006 /)), &
         plate( "KE", (/    6.4584,  39.9929,  2.3474 /)), &
         plate( "NI", (/ -174.4882,  -3.2883,  3.3136 /)), &
         plate( "TO", (/    4.4767,  25.8737,  8.9417 /)), &
         plate( "PM", (/ -113.9038,  31.3510,  0.3171 /)), &
         plate( "AS", (/  122.8665,  19.4251,  0.1239 /)), &
         plate( "AT", (/   26.6585,  40.1121,  1.2105 /)), &
         plate( "GP", (/   81.1806,   2.5287,  5.4868 /)), &
         plate( "EA", (/   67.5269,  24.9729, 11.3343 /)), &
         plate( "JZ", (/   70.7429,  34.2507, 22.3676 /)), &
         plate( "OK", (/  -92.2813,  30.3022,  0.2290 /)), &
         plate( "NB", (/  127.6370, -45.0406,  0.8563 /)), &
         plate( "SB", (/  -31.8883,   6.8767,  8.1107 /)), &
         plate( "MN", (/  150.2676,  -3.6699, 51.5690 /)), &
         plate( "NH", (/   -6.6018,   0.5684,  2.4688 /)), &
         plate( "BR", (/  142.0636, -63.7420,  0.4898 /)), &
         plate( "CL", (/   72.0525, -72.7849,  0.6066 /)), &
         plate( "MA", (/  137.8404,  11.0533,  1.3061 /)), &
         plate( "ND", (/ -122.6815,  17.7331,  0.1162 /)), &
         plate( "AP", (/  -83.9776,  -6.5763,  0.4881 /)), &
         plate( "BU", (/  -78.1008,  -6.1254,  2.2287 /)), &
         plate( "MS", (/  -56.0916,   2.1477,  3.5655 /)), &
         plate( "BS", (/  121.6413,  -1.4855,  2.4753 /)), &
         plate( "TI", (/  113.4976,  -4.4363,  1.8639 /)), &
         plate( "ON", (/  137.9182,  36.1163,  2.5391 /))  &
      /)                                                   &
   )                                                       &
/)
!-------------------------------------------------------------------------------
!  Plate name lookup table
      type(plate_names_list), parameter, private :: &
PLATE_NAMES = plate_names_list( & 
   (/ "AM", "AN", "AR", "AU", "CP", "CA", "CO", "EU", "IN", "JF", "LW", "MQ", &
      "NZ", "NA", "NB", "PA", "PS", "RI", "SW", "SC", "SM", "SA", "SU", "SR", &
      "YZ", "AS", "AP", "AT", "BR", "BS", "BH", "BU", "CL", "CR", "EA", "FT", &
      "GP", "JZ", "KE", "MN", "MO", "MA", "MS", "NH", "NI", "ND", "NB", "OK", &
      "ON", "PM", "SL", "SS", "SB", "TI", "TO", "WL" /), &
   (/ "Amur                          ", "Antarctica                    ",  &
		"Arabia                        ", "Australia                     ",  &
		"Capricorn                     ", "Caribbean                     ",  &
		"Cocos                         ", "Eurasia                       ",  &
		"India                         ", "Juan de Fuca                  ",  &
		"Lwandle                       ", "Macquarie                     ",  &
		"Nazca                         ", "North America                 ",  &
		"Nubia                         ", "Pacific                       ",  &
		"Philippine Sea                ", "Rivera                        ",  &
		"Sandwich                      ", "Scotia                        ",  &
		"Somalia                       ", "South America                 ",  &
		"Sunda                         ", "Sur                           ",  &
		"Yangtze                       ", "Aegean Sea                    ",  &
		"Altiplano                     ", "Anatolia                      ",  &
		"Balmoral Reef                 ", "Banda Sea                     ",  &
		"Birds Head                    ", "Burma                         ",  &
		"Caroline                      ", "Conway Reef                   ",  &
		"Easter                        ", "Futuna                        ",  &
		"Galapagos                     ", "Juan Fernandez                ",  &
		"Kermadec                      ", "Manus                         ",  &
		"Maoke                         ", "Mariana                       ",  &
		"Molucca Sea                   ", "New Hebrides                  ",  &
		"Niuafo'ou                     ", "North Andes                   ",  &
		"North Bismarck                ", "Okhotsk                       ",  &
		"Okinawa                       ", "Panama                        ",  &
		"Shetland                      ", "Solomon Sea                   ",  &
		"South Bismarck                ", "Timor                         ",  &
		"Tonga                         ", "Woodlark                      "   &
	/) &
)
      CONTAINS

!===============================================================================
   subroutine relative_plate_motion(lon,lat,P1,P2,model_name,v,az,rate)
!===============================================================================
!  Returns the spreading vector v(E,N) (mm/a) at a given point, when provided
!  with two plates.  The vector is for P1 wrt fixed P2.
!  Defaults to HS3-NUVEL1A
!  Can ask either for vector v, or for both az and rate.
   
      implicit none
      
      real(rs), intent(in) :: lon,lat
      real(rs)             :: v1(2),v2(2),rel(2)
      type(plate)  :: P1, P2
      character(len=*), optional, intent(in) :: model_name
      character(len=MODEL_NAME_LENGTH) :: model_name_in
      real(rs), optional, intent(out) :: v(2),az,rate
      
      model_name_in = 'HS3-NUVEL1A'
      
      if (present(model_name)) model_name_in = model_name
      
!  Get the rotation vectors for the plates in the desired model
      call get_plate(P1,trim(model_name_in))
      call get_plate(P2,trim(model_name_in))
      
!  Get the relative motion for the two plates: P1 wrt P2
      call absolute_plate_motion(lon,lat,P1,model_name=model_name_in,v=v1)
      call absolute_plate_motion(lon,lat,P2,model_name=model_name_in,v=v2)
      rel = v1 - v2
      
!  Return the desired values
      if (present(v)) v = rel
      if (present(az) .and. present(rate)) then
         az = mod(atan2(rel(1), rel(2))*180._rs/pi+3600._rs,360._rs)
         rate = sqrt(rel(1)**2 + rel(2)**2)
      else if ((present(az) .and. .not.present(rate)) .or. &
               (present(rate) .and. .not.present(az))) then
         write(lu_stderr,'(2a)') &
            'plate_motion: relative_plate_motion: Must ask for both az',&
            'and rate if requesting either.'
         stop
      endif
      
      return
   end subroutine relative_plate_motion
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine absolute_plate_motion(lon,lat,P,model_name,v,az,rate)
!===============================================================================
!  Returns the absolute plate motion for plate P at lon,lat
!  as a vector (E,N) (mm/a)

      implicit none
      
      real(rs) :: apm(2)
      real(rs), intent(in) :: lon,lat
      type(plate), intent(inout) :: P
      character(len=*), optional, intent(in) :: model_name
      character(len=MODEL_NAME_LENGTH) :: model_name_in
      real(rs), optional, intent(out) :: v(2),az,rate
      
      model_name_in = 'HS3-NUVEL1A'
      if (present(model_name)) model_name_in = model_name
      
!  Get the Euler pole for the plate
      call get_plate(P,trim(model_name))
      
!  Calculate the apm
      apm = calc_apm(lon,lat,P%pole(1),P%pole(2),P%pole(3))
!  Convert from km/Ma to mm/a
      apm = apm / 1000._rs
      
!  Return the desired values
      if (present(v)) v = apm
      if (present(az) .and. present(rate)) then
         az = mod(atan2(apm(1), apm(2))*180._rs/pi+3600._rs,360._rs)
         rate = sqrt(apm(1)**2 + apm(2)**2)
      else if ((present(az) .and. .not.present(rate)) .or. &
               (present(rate) .and. .not.present(az))) then
         write(lu_stderr,'(2a)') &
            'plate_motion: absolute_plate_motion: Must ask for both az',&
            'and rate if requesting either.'
         stop
      endif
  
      return
   end subroutine absolute_plate_motion
!-------------------------------------------------------------------------------

!===============================================================================
   function have_model(model_name)
!===============================================================================
!  Checks that the model is available in this module.

      implicit none
      
      character(len=*), intent(in) :: model_name
      logical :: have_model
      integer :: imodel
      
      have_model = .false.
      
!  Look for model in available models
      do imodel=1,size(MODEL)
         if (MODEL(imodel)%name == model_name) then
            have_model = .true.
            exit
         endif
      enddo
     
      return
   end function have_model
!-------------------------------------------------------------------------------


!===============================================================================
   subroutine get_plate(p,model_name)
!===============================================================================
!  Fill in the Euler pole values for the desired plate in the specified model,
!  whilst checking the model and the plates exist in the first place.

      implicit none
      
      type(plate), intent(inout) :: p
      character(len=*) :: model_name
      integer :: iplate,imodel
      
!  Check that we have this model
      if (.not.have_model(trim(model_name))) then
         write(lu_stderr,'(2a)') 'plate_motion: get_plate: Model not found: ',trim(model_name)
         stop
      endif
      
!  Find number of the model
      do imodel=1,size(MODEL)
         if (trim(MODEL(imodel)%name) == trim(model_name)) exit
      enddo
      
!  Check this plate exists 
      if (.not.any(MODEL(imodel)%p%name == p%name)) then
         write(lu_stderr,'(5a)') 'plate_motion: get_plate: Plate ',p%name,&
                                 ' not found in ',trim(MODEL(imodel)%name),' model.'
         stop
      endif
      
!  Get plate Euler pole and rotation rate
      do iplate=1,MODEL(imodel)%n_plates
         if (MODEL(imodel)%p(iplate)%name == p%name) then
            p%pole = MODEL(imodel)%p(iplate)%pole
            exit
         endif
      enddo
      
      return
   end subroutine get_plate
!-------------------------------------------------------------------------------

!===============================================================================
   function calc_apm(lon,lat,plon,plat,omega)
!===============================================================================
!  Calculate the motion of a point on a sphere given its lon,lat and and Euler pole
!  and rotation rate.  Output units are determined by input units: rotation
!  rates of deg/Ma will give m/Ma; deg/a will give m/a, etc.
!  Returns a vector pointing (E,N)

      implicit none
      
      real(rs), intent(in) :: lon,lat,plon,plat,omega
      real(rs) :: rlon,rlat,rplon,rplat,romega
      real(rs) :: calc_apm(2)
      real(rs) :: a(3), w(3), x(3), N(3), E(3)

!  Check for valid points on the sphere
      if (lat < -90._rs .or. lat > 90._rs .or. plat < -90._rs .or. plat > 90._rs ) then
         write(lu_stderr,'(a)') &
            'plate_motion: calc_apm: latitudes must be in the range -90 to 90 deg.'
         stop
      else if (lon < -180._rs .or. lon > 180._rs .or. plon < -180._rs .or. plon > 180._rs ) then
         write(lu_stderr,'(a)') &
            'plate_motion: calc_apm: longitudes must be in the range -180 to 180 deg.'
            stop
      endif
      
!  Convert to radians or radians/year
      rlon = lon * pi/180._rs  ;  rlat = lat * pi/180._rs
      rplon = plon*pi/180._rs  ;  rplat = plat*pi/180._rs
      romega = pi*omega/180._rs

!  Calculate rotation pole vector: latitude, not colatitude [/a or /Ma]
      w(1) = romega * cos(rplon) * cos(rplat)
      w(2) = romega * sin(rplon) * cos(rplat)
      w(3) = romega * sin(rplat)
      
!  Calculate location vector: latitude, not colatitude [m]
      x(1) = 1000._rs * R_EARTH * cos(rlon) * cos(rlat)
      x(2) = 1000._rs * R_EARTH * sin(rlon) * cos(rlat)
      x(3) = 1000._rs * R_EARTH * sin(rlat)
      
!  N-pointing unit vector [m]
      N(1) = -sin(rlat) * cos(rlon)
      N(2) = -sin(rlat) * sin(rlon)
      N(3) = cos(rlat)
      
!  E-pointing unit vector [m]
      E(1) = -sin(rlon)
      E(2) = cos(rlon)
      E(3) = 0._rs
      
!  Calculate motion vector = w (cross) x [m/a or m/Ma]
      a(1) = w(2)*x(3) - w(3)*x(2)
      a(2) = w(3)*x(1) - w(1)*x(3)
      a(3) = w(1)*x(2) - w(2)*x(1)
      
!  Find out components in plane normal to x, pointing N and 90 to it
      calc_apm(1) = dot_product(a,E)
      calc_apm(2) = dot_product(a,N)
      
      return
   end function calc_apm
!-------------------------------------------------------------------------------

!===============================================================================
   function get_plate_name(ab) result(long_name)
!===============================================================================
!  Gives you the long name of the plate, given the two-letter abbreviation.
   implicit none
   
   character(len=PLATE_NAME_LENGTH), intent(in) :: ab
   character(len=LONG_PLATE_NAME_LENGTH) :: long_name
   character(len=PLATE_NAME_LENGTH) :: ab_upper
   integer :: iplate
   
!  Convert to uppercase
      ab_upper = toupper(ab)
   
!  Go through plates in the PLATE_NAMES type
      do iplate = 1, size(PLATE_NAMES%ab)
         if (ab_upper == PLATE_NAMES%ab(iplate)) then
            long_name = PLATE_NAMES%name(iplate)
            exit
         endif
!  If we've reached here, no plate abbreviations match that given, so set long_name
!  to an indicative value
         long_name = NO_PLATE_NAME
      enddo
      
      if (long_name == NO_PLATE_NAME) &
         write(0,'(4a)') 'plate_motion: get_plate_name: ', &
                         'Warning: No plate with the abbreviation "',ab,'"'

   end function get_plate_name
!-------------------------------------------------------------------------------
    
!===============================================================================
   function toupper(string) result(upper)
!===============================================================================
!  Utility function to make a string uppercase
   implicit none
   
   character(len=*), intent(in) :: string
   character(len=len(string)) :: upper
   integer :: j
   
   do j = 1,len(string)
     if(string(j:j) >= "a" .and. string(j:j) <= "z") then
          upper(j:j) = achar(iachar(string(j:j)) - 32)
     else
          upper(j:j) = string(j:j)
     end if
   end do
   end function toupper
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine plate_list_all
!===============================================================================
!  Prints out information about which plates and models are available.
   implicit none
   
   integer :: i
   
   write(*,'(a)') '==============================================',&
                  'plate_motion: available plates and APM models:',&
                  '----------------------------------------------',&
                  'Plates:'
   do i = 1,size(PLATE_NAMES%ab)
      write(*,'("  ",a,2x,a)') PLATE_NAMES%ab(i),trim(PLATE_NAMES%name(i))
   enddo
   
   write(*,'(a)') '----------------------------------------------',&
                  'APM models:'
   do i = 1,N_AVAILABLE_MODELS
      write(*,'("  ",a)') trim(MODEL(i)%name)
   enddo
   
   write(*,'(a)') '==============================================='
   
   end subroutine plate_list_all
!-------------------------------------------------------------------------------

!_______________________________________________________________________________
end module plate_motion