//
//  smarties
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland. All rights reserved.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef smarties_Approximator_h
#define smarties_Approximator_h

#include "Builder.h"
#include "ThreadContext.h"
#include "../Utils/StatsTracker.h"
#include "../Math/Continuous_policy.h"
//#include "../Utils/SstreamUtilities.h"
//#include <iostream>
#include "../ReplayMemory/MemoryBuffer.h"

namespace smarties
{

class ParameterBlob;

struct Approximator
{
  //when this flag is true, specification of network properties is disabled:
  bool bCreatedNetwork = false;

  void setNumberOfAddedSamples(const Uint nSamples = 0);
  //specify type (and size) of auxiliary input
  void setAddedInput(const ADDED_INPUT type, Sint size = -1);
  // specify whether we are using target networks
  void setUseTargetNetworks(const Sint targetNetworkSampleID = -1,
                            const bool bTargetNetUsesTargetWeights = true);

  void setBlockGradsToPreprocessing();

  void buildFromSettings(const Uint outputSize) {
    buildFromSettings( std::vector<Uint>(1, outputSize) );
  }
  void buildFromSettings(const std::vector<Uint> outputSizes);
  void buildPreprocessing(const std::vector<Uint> outputSizes);
  Builder& getBuilder()
  {
    if (build) return * build.get();
    else {
      die("Requested unallocated network building entity");
      return * build.get();
    }
  }

  void initializeNetwork();

  Uint nOutputs() const {
    if (net == nullptr) return 0;
    else return net->getnOutputs();
  }
  Uint nLayers() const {
    if (net not_eq nullptr) return net->nLayers;
    else if (build not_eq nullptr) return build->layers.size();
    else return 0;
  }
  void setNgradSteps(const Uint iter) const { opt->nStep = iter; }
  void updateGradStats(const std::string& base, const Uint iter) const
  {
    gradStats->reduce_stats(base+"_"+name, iter);
  }

  Approximator(std::string name_, const HyperParameters&, const ExecutionInfo&,
               const MemoryBuffer* const replay_,
               const Approximator* const preprocessing_ = nullptr,
               const Approximator* const auxInputNet_ = nullptr);
  ~Approximator();

  void load(const MiniBatch& B, const Uint batchID, const Sint wghtID) const
  {
    const Uint thrID = omp_get_thread_num();
    // ensure we allocated enough workspaces:
    assert(contexts.size()>thrID && threadsPerBatch.size()>batchID);
    ThreadContext&C = * contexts[thrID].get();
    threadsPerBatch[batchID] = thrID;
    assert(C.endBackPropStep(0)<0 && "Previous backprop did not finish?");
    C.load(net, B, batchID, wghtID);
    //if(preprocessing) preprocessing->load(B, batchID, wghtID);
    //if(auxInputNet) auxInputNet->load(B, batchID, wghtID);
  }

  void load(const MiniBatch& B, const Agent& agent, const Sint wghtID=0) const
  {
    assert(agentsContexts.size() > agent.ID);
    AgentContext & C = * agentsContexts[agent.ID].get();
    C.load(net, B, agent, wghtID);
    //if(preprocessing) preprocessing->load(B, agent, wghtID);
    //if(auxInputNet) auxInputNet->load(B, agent, wghtID);
  }

  template< typename contextid_t, typename val_t>
  void setAddedInput(const std::vector<val_t>& addedInput,
                     const contextid_t& contextID,
                     const Uint t, Sint sampID = 0) const
  {
    getContext(contextID).setAddedInput(addedInput, t, sampID);
  }
  template< typename contextid_t>
  void setAddedInputType(const ADDED_INPUT& type,
                         const contextid_t& contextID,
                         const Uint t, Sint sampID = 0) const
  {
    getContext(contextID).setAddedInputType(type, sampID);
  }

  // forward: compute net output taking care also to gather additional required
  // inputs such as recurrent connections and auxiliary input networks.
  // It expects as input either the index over a previously loaded minibatch
  // or a previously loaded agent.
  template< typename contextid_t >
  Rvec forward(const contextid_t& contextID,
               const Uint t, Sint sampID=0, const bool overwrite=false) const
  {
    auto& C = getContext(contextID);
    if(sampID > (Sint) C.nAddedSamples) { sampID = 0; }
    if(overwrite)
       C.activation(t, sampID)->written = false;
    if(C.activation(t, sampID)->written)
      return C.activation(t, sampID)->getOutput();
    const Uint ind = C.mapTime2Ind(t);

    // Compute previous outputs if needed by recurrencies. LIMITATION:
    // What should we do for target net / additional samples?
    // next line assumes we want to use curr W and sample 0 for recurrencies
    // as well as using actions as added input instead of vectors/networks
    if(ind>0) assert(t>0);
    if(ind>0 && not C.activation(t-1, 0)->written) {
      const auto myInputType = C.addedInputType(0);
      if(myInputType not_eq NONE) C.addedInputType(0) = ACTION;
      forward(contextID, t-1, 0);
      if(myInputType not_eq NONE) C.addedInputType(0) = myInputType;
    }
    //if(ind>0 && not C.net(t, samp)->written) forward(C, t-1, samp);
    const Activation* const recur = ind>0? C.activation(t-1, 0) : nullptr;
    const Activation* const activation = C.activation(t, sampID);
    const Parameters* const W = opt->getWeights(C.usedWeightID(sampID));
    //////////////////////////////////////////////////////////////////////////////
    NNvec INP;
    if(preprocessing)
    {
      const Rvec preprocInp = preprocessing->forward(contextID, t, sampID);
      INP.insert(INP.end(), preprocInp.begin(), preprocInp.end());
    } else INP = C.getState(t);

    if(C.addedInputType(sampID) == NETWORK)
    {
      assert(auxInputNet not_eq nullptr);
      Rvec addedinp = auxInputNet->forward(contextID, t, sampID);
      assert( (Sint) addedinp.size() >= m_auxInputSize);
      addedinp.resize(m_auxInputSize);
      INP.insert(INP.end(), addedinp.begin(), addedinp.end());
    }
    else if(C.addedInputType(sampID) == ACTION)
    {
      const auto& addedinp = C.getAction(t);
      INP.insert(INP.end(), addedinp.begin(), addedinp.end());
    }
    else if(C.addedInputType(sampID) == VECTOR)
    {
      const auto& addedinp = C.addedInputVec(t, sampID);
      INP.insert(INP.end(), addedinp.begin(), addedinp.end());
    }
    assert(INP.size() == net->getnInputs());
    ////////////////////////////////////////////////////////////////////////////
    return net->forward(INP, recur, activation, W);
  }

  // forward target network
  template< typename contextid_t >
  Rvec forward_tgt(const contextid_t& contextID,
                   const Uint t, const bool overwrite=false) const
  {
    return forward(contextID, t, -1, overwrite);
  }

  // run network for agent's recent step
  Rvec forward(const Agent& agent, const bool overwrite = false) const
  {
    const auto& C = getContext(agent);
    return forward(agent, C.episode()->nsteps() - 1, 0, overwrite);
  }

  void setGradient(const Rvec& gradient,
                   const Uint  batchID,
                   const Uint  t, Sint sampID = 0) const
  {
    ThreadContext& C = getContext(batchID);
    if(sampID > (Sint) C.nAddedSamples) { sampID = 0; }
    //for(Uint i=0; i<grad.size(); ++i) grad[i] *= PERW;
    gradStats->track_vector(gradient, C.threadID);
    const Sint ind = C.mapTime2Ind(t);
    //ind+1 because we use c-style for loops in other places:
    C.endBackPropStep(sampID) = std::max(C.endBackPropStep(sampID), ind+1);
    assert( C.activation(t, sampID)->written );
    if(ESpopSize > 1) debugL("Skipping backward because we use ES.");
    else C.activation(t, sampID)->addOutputDelta(gradient);
  }

  Real& ESloss(const Uint ESweightID = 0) { return losses[ESweightID]; }
  Rvec oneStepBackProp(const Rvec& gradient,
                       const Uint  batchID,
                       const Uint  t, Sint sampID) const
  {
    assert(auxInputNet && "improperly set up the aux input net");
    assert(auxInputAttachLayer >= 0 && "improperly set up the aux input net");
    if(ESpopSize > 1) {
      debugL("Skipping backprop because we use ES optimizers.");
      return Rvec(m_auxInputSize, 0);
    }
    ThreadContext& C = getContext(batchID);
    if(sampID > (Sint) C.nAddedSamples) { sampID = 0; }
    const MDPdescriptor & MDP = replay->MDP;
    const Parameters* const W = opt->getWeights(C.usedWeightID(sampID));
    const Uint inputSize = preprocessing? preprocessing->nOutputs()
                         : (1+MDP.nAppendedObs) * MDP.dimStateObserved;
    Activation* const A = C.activation(t, sampID);
    assert(C.activation(t, sampID)->written == true);
    const Rvec ret = net->backPropToLayer(gradient, auxInputAttachLayer, A, W);
    //C.endBackPropStep(sampID) = -1; //to stop additional backprops
    //printf("%f\n", ret[inputSize]);
    if(auxInputAttachLayer>0) return ret;
    else return Rvec(& ret[inputSize], & ret[inputSize + m_auxInputSize]);
  }

  Rvec getStepBackProp(const Uint batchID, const Uint t, Sint sampID) const
  {
    assert(auxInputNet && "improperly set up the aux input net");
    assert(auxInputAttachLayer >= 0 && "improperly set up the aux input net");
    if(ESpopSize > 1) {
      debugL("Skipping backprop because we use ES optimizers.");
      return Rvec(m_auxInputSize, 0);
    }
    ThreadContext& C = getContext(batchID);
    if(sampID > (Sint) C.nAddedSamples) { sampID = 0; }
    const MDPdescriptor & MDP = replay->MDP;
    const Uint inputSize = preprocessing? preprocessing->nOutputs()
                         : (1+MDP.nAppendedObs) * MDP.dimStateObserved;
    assert(C.activation(t, sampID)->written == true);
    Rvec ret = C.activation(t, sampID)->getInputGradient(auxInputAttachLayer);
    //C.endBackPropStep(sampID) = -1; //to stop additional backprops
    //printf("%f\n", ret[inputSize]);
    if(auxInputAttachLayer>0) return ret;
    else return Rvec(& ret[inputSize], & ret[inputSize + m_auxInputSize]);
  }

  void backProp(const Uint batchID) const
  {
    ThreadContext& C = getContext(batchID);

    if(ESpopSize > 1)
    {
      debugL("Skipping gradient because we use ES (derivative-free) optimizers.");
    }
    else
    {
      const auto& activations = C.activations;
      const auto& timeSeriesBase = activations[0];

      //loop over all the network samples, each may need different BPTT window
      for(Uint j = 0; j < activations.size(); ++j)
      {
        const Uint samp = activations.size() - 1 - j;
        const Sint last_error = C.endBackPropStep(samp);
        if(last_error < 0) continue;

        const auto& timeSeriesSamp = activations[samp];
        for (Sint i=0; i<last_error; ++i) assert(timeSeriesSamp[i]->written);

        const Parameters* const W = opt->getWeights(C.usedWeightID(samp));
        net->backProp(timeSeriesSamp, timeSeriesBase, last_error,
                      C.partialGradient.get(), W);

        if(preprocessing and not m_blockInpGrad)
        {
          for(Sint k=0; k<last_error; ++k)
          {
            const Uint t = C.mapInd2Time(k);
            // assume that preprocessing is layer 0:
            Rvec inputGrad = C.activation(t, samp)->getInputGradient(0);
            // we might have added inputs, therefore trim those:
            inputGrad.resize(preprocessing->nOutputs());
            preprocessing->setGradient(inputGrad, batchID, t, samp);
          }
        }
        C.endBackPropStep(samp) = -1; //to stop additional backprops
      }
    }

    nAddedGradients++;
  }

  void prepareUpdate()
  {
    #ifndef NDEBUG
    for(const auto& C : contexts)
      for(const Sint todoBackProp : C->lastGradTstep)
        assert(todoBackProp<0 && "arrived into prepareUpdate() before doing backprop on all workspaces");
    #endif

    if(nAddedGradients==0) die("No-gradient update. Revise hyperparameters.");

    if(preprocessing and not m_blockInpGrad)
      for(Uint i=0; i<ESpopSize; ++i) preprocessing->losses[i] += losses[i];

    opt->prepare_update(losses);
    losses = Rvec(ESpopSize, 0);
    reducedGradients = 1;
    nAddedGradients = 0;
  }

  bool ready2ApplyUpdate()
  {
    if(reducedGradients == 0) return true;
    else return opt->ready2UpdateWeights();
  }

  void applyUpdate()
  {
    if(reducedGradients == 0) return;

    opt->apply_update();
    reducedGradients = 0;
  }

  void gatherParameters(ParameterBlob& params) const;

  void getHeaders(std::ostringstream& buff) const;
  void getMetrics(std::ostringstream& buff) const;
  void save(const std::string base, const bool bBackup);
  void restart(const std::string base);
  void rename(std::string newname) { name = newname; }

  Optimizer * getOptimizerPtr() {
    return opt.get();
  }
  Network * getNetworkPtr() {
    return net.get();
  }

private:
  const HyperParameters & settings;
  const ExecutionInfo & distrib;
  std::string name;
  const Uint   nAgents =  distrib.nAgents,    nThreads =  distrib.nThreads;
  const Uint ESpopSize = settings.ESpopSize, batchSize = settings.batchSize;
  const MemoryBuffer* const replay;
  const Approximator* const preprocessing;
  const Approximator* const auxInputNet;
  const ActionInfo& aI = replay->aI;
  Sint auxInputAttachLayer = -1;
  Sint m_auxInputSize = -1;
  Uint m_numberOfAddedSamples = 0;
  bool m_UseTargetNetwork = false;
  bool m_bTargetNetUsesTargetWeights = true;
  Sint m_targetNetworkSampleID = -1;

  // Whether to backprop gradients in the input network.
  // Some papers do not propagate policy gradients towards encoding layers
  bool m_blockInpGrad = false;

  std::shared_ptr<Network> net;
  std::shared_ptr<Optimizer> opt;
  std::unique_ptr<Builder> build;

  mutable std::vector<Uint> threadsPerBatch = std::vector<Uint>(batchSize, -1);
  std::vector<std::unique_ptr<ThreadContext>> contexts;
  std::vector<std::unique_ptr< AgentContext>> agentsContexts;
  StatsTracker* gradStats = nullptr;

  // For CMAES based optimization. Keeps track of total loss associate with
  // Each weight vector sample:
  mutable Rvec losses = Rvec(ESpopSize, 0);

  ThreadContext& getContext(const Uint batchID) const
  {
    assert(threadsPerBatch.size() > batchID);
    return * contexts[ threadsPerBatch[batchID] ].get();
  }
  AgentContext&  getContext(const Agent& agent) const
  {
    assert(agentsContexts.size() > agent.ID && agentsContexts[agent.ID]);
    return * agentsContexts[agent.ID].get();
  }

public:
  mutable std::atomic<Uint> nAddedGradients{0};
  Uint reducedGradients=0;
};

} // end namespace smarties
#endif // smarties_Approximator_h
