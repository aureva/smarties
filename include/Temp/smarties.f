!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c    This program is developed for parallel            c
!     c    computation of scalar transport in ABL            c
!     c                                                      c
!     c					    July, 2002       c
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      Program  smarties
 
      implicit none

      interface
      subroutine smarties_sendInitState( 
     +   ptr2comm, state, state_dim, agentID) bind(c, 
     +   name='smarties_sendInitState')
      
      use, intrinsic :: iso_c_binding
      
      implicit none
      
      type(c_ptr),    intent(in), value :: ptr2comm
      type(c_ptr),    intent(in), value :: state
      integer(c_int), intent(in), value :: state_dim
      integer(c_int), intent(in), optional :: agentID
      
      end subroutine smarties_sendInitState
      end interface

      interface
      subroutine smarties_sendState( 
     +   ptr2comm, state, state_dim, reward, agentID) 
     +   bind(c, name='smarties_sendState')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr),    intent(in), value :: ptr2comm
      type(c_ptr),    intent(in), value :: state
      integer(c_int), intent(in), value :: state_dim
      real(c_double), intent(in), value :: reward
      integer(c_int), intent(in), optional :: agentID
      end subroutine smarties_sendState
      end interface
            
      interface
      subroutine smarties_sendLastState( 
     +   ptr2comm, state, state_dim, reward, agentID) 
     +   bind(c, name='smarties_sendLastState')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr),    intent(in), value :: ptr2comm
      type(c_ptr),    intent(in), value :: state
      integer(c_int), intent(in), value :: state_dim
      real(c_double), intent(in), value :: reward
      integer(c_int), intent(in), optional :: agentID
      end subroutine smarties_sendLastState
      end interface
 
      interface
      subroutine smarties_recvAction( 
     +   ptr2comm, action, action_dim, agentID) 
     +   bind(c, name='smarties_recvAction')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr),    intent(in), value :: ptr2comm
      type(c_ptr),    intent(in), value :: action
      integer(c_int), intent(in), value :: action_dim
      integer(c_int), intent(in), optional :: agentID
      end subroutine smarties_recvAction
      end interface

      interface
      subroutine smarties_setNumAgents( 
     +   ptr2comm, num_agents) 
     +   bind(c, name='smarties_setNumAgents')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr),    intent(in), value :: ptr2comm
      integer(c_int), intent(in), value :: num_agents
      end subroutine smarties_setNumAgents
      end interface

      interface
      subroutine smarties_setActionScales( 
     +   ptr2comm, upper_act_bound, lower_act_bound, 
     +   bounded, action_dim, agentID) 
     +   bind(c, name='smarties_setActionScales')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr),     intent(in), value :: ptr2comm
      type(c_ptr),     intent(in), value :: upper_act_bound
      type(c_ptr),     intent(in), value :: lower_act_bound
      logical(c_bool), intent(in), value :: bounded
      integer(c_int),  intent(in), value :: action_dim
      integer(c_int),  intent(in), optional :: agentID
      end subroutine smarties_setActionScales
      end interface

      interface
      subroutine smarties_setStateActionDims( 
     +   ptr2comm, state_dim, action_dim, agentID) 
     +   bind(c, name='smarties_setStateActionDims')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr),    intent(in), value :: ptr2comm
      integer(c_int), intent(in), value :: state_dim
      integer(c_int), intent(in), value :: action_dim
      integer(c_int), intent(in), optional :: agentID
      end subroutine smarties_setStateActionDims
      end interface

      interface
      subroutine smarties_setStateObservable( 
     +   ptr2comm, b_observable, state_dim, agentID) 
     +   bind(c, name='smarties_setStateObservable')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr),    intent(in), value :: ptr2comm
      type(c_ptr),    intent(in), value :: b_observable
      integer(c_int), intent(in), value :: state_dim
      integer(c_int), intent(in), optional :: agentID
      end subroutine smarties_setStateObservable
      end interface
 
      END program smarties

!=================================================================
