!======================================================================================================================!
!                                                                                                                      !
!     Program:      CANACC                                                                                             !
!                   Canopy physics model for ACCESS                                                                    !
!                                                                                                                      !
!     Version:      2.0.0                                                                                              !
!                                                                                                                      !
!     Contact:      Rick D. Saylor, PhD                                                                                !
!                   Physical Scientist                                                                                 !
!                   U. S. Department of Commerce                                                                       !
!                   National Oceanic and Atmospheric Administration                                                    !
!                   Air Resources Laboratory                                                                           !
!                   Atmospheric Turbulence and Diffusion Division                                                      !
!                   456 S. Illinois Ave                                                                                !
!                   Oak Ridge, TN 37830                                                                                !
!                   email: Rick.Saylor@noaa.gov                                                                        !
!                                                                                                                      !
!======================================================================================================================!
!**********************************************************************************************************************!
!************************************************ LEGAL NOTICE ********************************************************!
!**********************************************************************************************************************!
!                                                                                                                      ! 
!     This computer software is in the public domain and is not protected by U. S. copyright law as it was created     !
!     by a U. S. Federal Government employee in the normal course of their official duties. However, the author (RDS)  ! 
!     makes the following notifications and requests:                                                                  ! 
!                                                                                                                      ! 
!     (1) Please report any found software bugs, problems or unusual results to "Rick.Saylor@noaa.gov".                !
!     (2) Do not redistribute the CANACC code, peripheral tools, or its direct derivatives without notifying the       !
!         author at "Rick.Saylor@noaa.gov".                                                                            ! 
!     (3) Neither the U. S. Government nor the author of CANACC is obligated to maintain, provide updates to, fix      !
!         defects in or continue development of the CANACC code or peripheral tools.                                   !
!     (4) The CANACC code and peripheral tools (the "SOFTWARE") are subject to the LIMITATION OF LIABILITY as stated   !
!         below.                                                                                                       ! 
!     (5) The CANACC code and peripheral tools are only intended for scientific research or educational purposes and   !
!         should not be used for or in any commercial application.                                                     !
!     (6) Any works derived from the CANACC source code or peripheral tools are the sole responsibility of the         !
!         modifying developer and not of the U. S. Government or the original author of CANACC.                        ! 
!     (7) Any modifications, extensions or additions to the CANACC code and peripheral tools must strictly follow      !
!         the coding style and conventions described in the "CANACC User's Guide and Documentation" distributed        !
!         with CANACC.                                                                                                 ! 
!                                                                                                                      ! 
!     LIMITATION OF LIABILITY: THE UNITED STATES GOVERNMENT DOES NOT WARRANT THE ACCURACY OR COMPLETENESS OF THIS      !
!     SOFTWARE, AND IS NOT RESPONSIBLE FOR ERRORS AND/OR OMISSIONS, AND IS NOT LIABLE FOR ANY DIRECT, INDIRECT,        ! 
!     INCIDENTAL, CONSEQUENTIAL, SPECIAL OR EXEMPLARY DAMAGES OR LOST PROFIT RESULTING FROM ANY USE OR MISUSE OF       ! 
!     THIS SOFTWARE OR THE INFORMATION GENERATED BY IT. THIS SOFTWARE IS BEING DISTRIBUTED "AS IS" AND THE U.S.        !
!     GOVERNMENT DOES NOT MAKE ANY WARRANTY CLAIMS, EITHER EXPRESSED OR IMPLIED, WITH RESPECT TO ITS QUALITY,          !
!     ACCURACY, COMPLETENESS, PERFORMANCE, MERCHANTABILITY OR FITNESS FOR ANY INTENDED PURPOSE.                        ! 
!                                                                                                                      ! 
!**********************************************************************************************************************!
!======================================================================================================================!
program CANACC
   use GlobalData
   use Initialize
   use EnvironData
   use CanopyPhysics
   use Output
   use Utils
   implicit none
   real(kind=4) :: t0, tf

   call cpu_time(t0)

   call InitializeModel()

   do 
      ! read environmental data
      call ReadEnvironData()

      ! partition measured solar radiation and PPFD into direct and diffuse components
      call PartitionRAD()

      ! calculate profiles through the entire domain
      call CalcRadProfiles()

      ! calculate canopy physics, including temperature profiles, photosynthetic assimilation rates, and
      ! stomatal conductances
      call CalcCanopyPhysics()

      ! calculate sun/shade weighted canopy profiles
      call CalcWeightedProfiles()
      
      ! output summary results and store everything
      call OutputResult()

      nt=nt+1
      t = t + dtout

      if (t >= tend) exit
   end do 

   call PrintFinaltoFile()

   call CleanUp()

   call cpu_time(tf)

   call PrintCPUtime((tf-t0))

end program CANACC
