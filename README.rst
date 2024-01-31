# smarties
#coupled with LESGO

Reinforcement learning wall-model (RLWM) for large eddy simulation
smarties library coupled with LESGO

-> You should first download and install smarties from https://github.com/cselab/smarties

-> The /apps folders contains the routines to run the model, for example LESGO with coupling routines (smarties_stat.f90 and app_main.f90). In particular, the functions send_recv_state_action is the main important since it enables to communicate with agents and get the action. 
CAREFUL : use runfile_training with care since it erases the existing weights and biases in the /runs folder if it exists.
Use runfile_exec to use the RLWM without retraining.
 Details can be found in :
Vadrot, A., Yang, X. I., Bae, H. J., & Abkar, M. (2023). Log-law recovery through reinforcement-learning wall model for large eddy simulation. Physics of Fluids, 35(5).

-> You can execute the code in two ways: either train your own model or run the existing model. In case you train your own model, it will create a new folder in /runs and print the output of LESGO for all the simulations that are used for training simulation_000_....  When testing it will use the weights and biases of RL network (agent_00_net***) contained in the corresponding /runs folder. 
CAREFUL: use runfile_training with care since it erases the existing weights and biases in the /runs folder if they exist.  

-> To cite this repository, reference the paper:
@article{vadrot2023log,
  title={Log-law recovery through reinforcement-learning wall model for large eddy simulation},
  author={Vadrot, Aur{\'e}lien and Yang, Xiang IA and Bae, H Jane and Abkar, Mahdi},
  journal={Physics of Fluids},
  volume={35},
  number={5},
  year={2023},
  publisher={AIP Publishing}
}

