//
//  smarties
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland. All rights reserved.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "Communicator.h"
#include "Utils/SocketsLib.h"
#include "Core/Worker.h"

namespace smarties
{

void Communicator::setStateActionDims(const int dimState,
                                         const int dimAct,
                                         const int agentID)
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  if( (size_t) agentID >= ENV.descriptors.size())
    die("Attempted to write to uninitialized MDPdescriptor");

  ENV.descriptors[agentID]->dimState = dimState;
  ENV.descriptors[agentID]->dimAction = dimAct;
}

void Communicator::setActionScales(const std::vector<double> uppr,
                                     const std::vector<double> lowr,
                                     const bool bound,
                                     const int agentID)
{
  setActionScales(uppr,lowr, std::vector<bool>(uppr.size(),bound), agentID);
}
void Communicator::setActionScales(const std::vector<double> upper,
                                   const std::vector<double> lower,
                                   const std::vector<bool>   bound,
                                   const int agentID)
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  if(agentID >= (int) ENV.descriptors.size())
    die("Attempted to write to uninitialized MDPdescriptor");
  if(upper.size() not_eq ENV.descriptors[agentID]->dimAction or
     lower.size() not_eq ENV.descriptors[agentID]->dimAction or
     bound.size() not_eq ENV.descriptors[agentID]->dimAction )
    die("size mismatch");
  if(ENV.descriptors[agentID]->bDiscreteActions())
    die("either continuous or discrete actions");

  ENV.descriptors[agentID]->upperActionValue =
                Rvec(upper.begin(), upper.end());
  ENV.descriptors[agentID]->lowerActionValue =
                Rvec(lower.begin(), lower.end());
  ENV.descriptors[agentID]->bActionSpaceBounded =
    std::vector<bool>(bound.begin(), bound.end());
}

void Communicator::setActionOptions(const int options,
                                      const int agentID)
{
  setActionOptions(std::vector<int>(1, options), agentID);
}

void Communicator::setActionOptions(const std::vector<int> options,
                                      const int agentID)
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  if(agentID >= (int) ENV.descriptors.size())
    die("Attempted to write to uninitialized MDPdescriptor");
  if(options.size() not_eq ENV.descriptors[agentID]->dimAction)
    die("size mismatch");

  ENV.descriptors[agentID]->discreteActionValues =
    std::vector<Uint>(options.begin(), options.end());
}

void Communicator::setStateObservable(const std::vector<bool> observable,
                                        const int agentID)
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  if(agentID >= (int) ENV.descriptors.size())
    die("Attempted to write to uninitialized MDPdescriptor");
  if(observable.size() not_eq ENV.descriptors[agentID]->dimState)
    die("size mismatch");

  ENV.descriptors[agentID]->bStateVarObserved =
    std::vector<bool>(observable.begin(), observable.end());
}

void Communicator::setStateScales(const std::vector<double> upper,
                                    const std::vector<double> lower,
                                    const int agentID)
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state.");
    return;
  }
  if(agentID >= (int) ENV.descriptors.size())
    die("Attempted to write to uninitialized MDPdescriptor");
  const Uint dimS = ENV.descriptors[agentID]->dimState;
  if(upper.size() not_eq dimS or lower.size() not_eq dimS )
    die("size mismatch");

  // For consistency with action space we ask user for a rough box of state vars
  // but in reality we scale with mean and stdev computed during training.
  // This function serves only as an optional initialization for statistiscs.
  NNvec meanState(dimS), diffState(dimS);
  for (Uint i=0; i<dimS; ++i) {
    meanState[i] = (upper[i]+lower[i])/2;
    diffState[i] = std::fabs(upper[i]-lower[i]);
  }
  ENV.descriptors[agentID]->stateMean   = meanState;
  ENV.descriptors[agentID]->stateStdDev = diffState;
}

void Communicator::setIsPartiallyObservable(const int agentID)
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  if(agentID >= (int) ENV.descriptors.size())
    die("Attempted to write to uninitialized MDPdescriptor");

  ENV.descriptors[agentID]->isPartiallyObservable = true;
}

void Communicator::setPreprocessingConv2d(
  const int input_width, const int input_height, const int input_features,
  const int kernels_num, const int filters_size, const int stride,
  const int agentID)
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  if(agentID >= (int) ENV.descriptors.size())
    die("Attempted to write to uninitialized MDPdescriptor");

  // can be made to be more powerful (different sizes in x/y, padding, etc)
  Conv2D_Descriptor descr;
  descr.inpFeatures = input_features;
  descr.inpY        = input_height;
  descr.inpX        = input_width;
  descr.outFeatures = kernels_num;
  descr.filterx     = filters_size;
  descr.filtery     = filters_size;
  descr.stridex     = stride;
  descr.stridey     = stride;
  descr.paddinx     = 0;
  descr.paddiny     = 0;
  descr.outY   = (descr.inpY -descr.filterx +2*descr.paddinx)/descr.stridex + 1;
  descr.outX   = (descr.inpX -descr.filtery +2*descr.paddiny)/descr.stridey + 1;
  ENV.descriptors[agentID]->conv2dDescriptors.push_back(descr);
}

void Communicator::setNumAppendedPastObservations(
  const int n_appended, const int agentID)
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  if(agentID >= (int) ENV.descriptors.size())
    die("Attempted to write to uninitialized MDPdescriptor");

  ENV.descriptors[agentID]->nAppendedObs = n_appended;
}

void Communicator::setNumAgents(int _nAgents)
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  assert(_nAgents > 0);

  ENV.nAgentsPerEnvironment = _nAgents;
}

void Communicator::envHasDistributedAgents()
{
  /*
  if(comm_inside_app == MPI_COMM_NULL) {
    printf("ABORTING: Distributed agents has no effect on single-process "
    " applications. It means that each simulation rank holds different agents.");
    fflush(0); abort();
    bEnvDistributedAgents = false;
    return;
  }
  */
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  if(ENV.bAgentsHaveSeparateMDPdescriptors)
    die("ABORTING: Smarties supports either distributed agents (ie each "
    "worker holds some of the agents) or each agent defining a different MDP "
    "(state/act spaces).");

  bEnvDistributedAgents =  true;
}

void Communicator::agentsDefineDifferentMDP()
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  if(bEnvDistributedAgents) {
    printf("ABORTING: Smarties supports either distributed agents (ie each "
    "worker holds some of the agents) or each agent defining a different MDP "
    "(state/act spaces)."); fflush(0); abort();
  }

  ENV.initDescriptors(true);
}

void Communicator::disableDataTrackingForAgents(int agentStart, int agentEnd)
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }

  ENV.bTrainFromAgentData.resize(ENV.nAgentsPerEnvironment, 1);
  for(int i=agentStart; i<agentEnd; ++i)
    ENV.bTrainFromAgentData[i] = 0;
}

void Communicator::agentsShareExplorationNoise(const int agentID)
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  if(bEnvDistributedAgents) {
    printf("ABORTING: Shared exploration noise is not yet supported for "
      "distributed agents (noise is not communicated with MPI to other ranks "
      "running the simulation)."); fflush(0); abort();
  }
  if(agentID >= (int) ENV.descriptors.size())
    die("Attempted to write to uninitialized MDPdescriptor");
  ENV.descriptors[agentID]->bAgentsShareNoise = true;
}

void Communicator::finalizeProblemDescription()
{
  if(ENV.bFinalized) {
    warn("Cannot edit env description after having sent first state."); return;
  }
  synchronizeEnvironments();
}

void Communicator::_sendState(const int agentID, const episodeStatus status,
    const std::vector<double>& state, const double reward)
{
  //printf("Communicator.cpp sendState 1");
  //fflush(stdout);
  //if( not ENV.bFinalized ) printf("sendState 1 bis"); // race condition
  //if( not ENV.bFinalized ) fflush(stdout); // race condition
  if( not ENV.bFinalized ) synchronizeEnvironments(); // race condition
  if(bTrainIsOver)
    die("App recvd end-of-training signal but did not abort on it's own.");

  //printf("sendState 2");
  //fflush(stdout);

  //const auto& MDP = ENV.getDescriptor(agentID);
  assert(agentID>=0 && (Uint) agentID < agents.size());
  assert(agents[agentID]->localID == (unsigned) agentID);
  assert(agents[agentID]->ID == (unsigned) agentID);
  agents[agentID]->update(status, state, reward);

  //printf("sendState 3");
  //fflush(stdout);

  //#ifndef NDEBUG
  if (agents[agentID]->stateIsInvalid())
    die("Environment gave a nan or inf state or reward.");
  //#endif

  if(SOCK.server == -1)
  {
    assert(worker not_eq nullptr);
    worker->stepWorkerToMaster( * agents[agentID].get() );
  }
  else
  {
    agents[agentID]->packStateMsg(BUFF[agentID]->dataStateBuf);
    SOCKET_Bsend(BUFF[agentID]->dataStateBuf,
                 BUFF[agentID]->sizeStateMsg,
                 SOCK.server);
    SOCKET_Brecv(BUFF[agentID]->dataActionBuf,
                 BUFF[agentID]->sizeActionMsg,
                 SOCK.server);
    agents[agentID]->unpackActionMsg(BUFF[agentID]->dataActionBuf);
  }

  if(status >= LAST) {
    agents[agentID]->learnerAvgCumulativeReward = agents[agentID]->action[0];
  }
  // we cannot control application. if we received a termination signal we abort
  if(agents[agentID]->learnStatus == KILL) {
    printf("App recvd end-of-training signal.\n");
    bTrainIsOver = true;
  }
}

const std::vector<double> Communicator::recvAction(const int agentID) const
{
  assert( agents[agentID]->agentStatus < LAST && "Application read action for "
    "a terminal state or truncated episode. Undefined behavior.");
  return agents[agentID]->getAction();
}

int Communicator::recvDiscreteAction(const int agentID) const
{
  assert( agents[agentID]->agentStatus < LAST && "Application read action for "
    "a terminal state or truncated episode. Undefined behavior.");
  return (int) agents[agentID]->getDiscreteAction();
}

void Communicator::synchronizeEnvironments()
{
  if ( ENV.bFinalized ) return;

  if(SOCK.server == -1)
  {
    //printf("synchronize env A");
    //fflush(stdout);
    assert(worker not_eq nullptr);

    worker->synchronizeEnvironments();

  }
  else
  {
    initOneCommunicationBuffer();
    const auto sendBufferFunc = [&](void* buffer, size_t size) {
      SOCKET_Bsend(buffer, size, SOCK.server);
    };
    ENV.synchronizeEnvironments(sendBufferFunc);

    // allocate rest of communication buffers:
    for(size_t i=1; i<agents.size(); ++i) initOneCommunicationBuffer();
  }
  assert(BUFF.size() > 0);
}

void Communicator::initOneCommunicationBuffer()
{
  Uint maxDimState  = 0, maxDimAction = 0;
  assert(ENV.descriptors.size() > 0);
  for(size_t i=0; i<ENV.descriptors.size(); ++i)
  {
    maxDimState  = std::max(maxDimState,  ENV.descriptors[i]->dimState );
    maxDimAction = std::max(maxDimAction, ENV.descriptors[i]->dimAction);
  }
  assert(ENV.nAgentsPerEnvironment>0);
  assert(maxDimAction>0); // state can be 0-D
  BUFF.emplace_back(std::make_unique<COMM_buffer>(maxDimState, maxDimAction) );
}

std::mt19937& Communicator::getPRNG() {
  return gen;
}
Real Communicator::getUniformRandom(const Real begin, const Real end)
{
  std::uniform_real_distribution<Real> distribution(begin, end);
  return distribution(gen);
}
Real Communicator::getNormalRandom(const Real mean, const Real stdev)
{
  std::normal_distribution<Real> distribution(mean, stdev);
  return distribution(gen);
}

bool Communicator::isTraining() const {
  return bTrain;
}
bool Communicator::terminateTraining() const {
  return bTrainIsOver;
}

unsigned Communicator::getLearnersGradStepsNum(const int agentID)
{
  return agents[agentID]->learnerGradStepID;
}
unsigned Communicator::getLearnersTrainingTimeStepsNum(const int agentID)
{
  return agents[agentID]->learnerTimeStepID;
}
double Communicator::getLearnersAvgCumulativeReward(const int agentID)
{
  return agents[agentID]->learnerAvgCumulativeReward;
}

Communicator::Communicator(Worker*const W, std::mt19937&G, bool isTraining) :
gen(G()), bTrain(isTraining), worker(W) {}

}
