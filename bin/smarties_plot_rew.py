#!/usr/bin/env python3
#
#  smarties
#  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland. All rights reserved.
#  Distributed under the terms of the MIT license.
#
#  Created by Guido Novati (novatig@ethz.ch).
#
#
import argparse, sys, time, numpy as np, os, matplotlib.pyplot as plt, glob
from six.moves import input

colors = ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a', '#b15928', \
          '#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6', '#ffff99']

def getWorkerRankToPlot(parsed, agentid, PATH):
    if(parsed.workerRank < 0):
        for rank in range(2): # defaults to either 0 or 1
            FILE = "%s/agent_%02d_rank_%03d_cumulative_rewards.dat" \
                   % (PATH, agentid, rank)
            if os.path.isfile(FILE): return rank
        return -1 # nothing found
    else: return parsed.workerRank # return user's choice

def plotReward(args, ax, ensambleList, agentid, rank, colorid):
    L = parsed.averagingWindow
    allDataX, allDataY = None, None

    for PATH in ensambleList:
        FILE = "%s/agent_%02d_rank_%03d_cumulative_rewards.dat" % (PATH,agentid,rank)
        assert os.path.isfile(FILE)
        if args.bAsk and input("Display file %s? (y/n) " % (FILE)) == 'n' :
           continue

        DATA = np.fromfile(FILE, sep=' ')
        DATA = DATA.reshape(DATA.size//5, 5) # t_step grad_step worker seqlen R
        N = (DATA.shape[0] // L)*L
        span = DATA.shape[0] - N + np.arange(0, N)
        X = DATA[span, 1] - DATA[0,1]
        #X = np.arange(0, N)
        Y = DATA[span, 4]
        newSize0 = X.size//L
        if newSize0 == 0 : continue
        X = X.reshape(newSize0, L)
        Y = Y.reshape(newSize0, L)
        if allDataX is None:
            allDataX, allDataY = X, Y
        else:
            oldSize0, oldSize1 = allDataX.shape[0], allDataX.shape[1]
            size0 = min(oldSize0, X.shape[0])
            size1 = oldSize1 + L
            newAllDataX = np.zeros((size0, size1))
            newAllDataY = np.zeros((size0, size1))
            newAllDataX[:size0, :oldSize1] = allDataX[:size0, :]
            newAllDataY[:size0, :oldSize1] = allDataY[:size0, :]
            newAllDataX[:size0, oldSize1:] = X[:size0, :]
            newAllDataY[:size0, oldSize1:] = Y[:size0, :]
            allDataX = newAllDataX
            allDataY = newAllDataY


    X = allDataX.mean(1)
    Yb = np.percentile(allDataY, 20, axis=1)
    Yt = np.percentile(allDataY, 80, axis=1)
    #Ym = np.percentile(Y, 50, axis=1)
    Ym = np.mean(allDataY, axis=1)
    legend = "%s agent %d" % (PATH, agentid)
    fill  = ax.fill_between(X,Yb,Yt, facecolor=colors[colorid], alpha=0.5)
    line, = ax.plot(X,Ym, color=colors[colorid], label=legend)
    return fill, line

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Smarties reward plotting tool.")
    parser.add_argument('folders', nargs='+',
      help="List of run directories from which to plot the cumulative rewards")
    parser.add_argument('--averagingWindow', type=int, default=50,
      help="Number of episodes used used to compute cumulative " \
           "rewards averages and standard deviations.")
    parser.add_argument('--workerRank', type=int, default=-1,
      help="Will plot the cumulative rewards obtained by this MPI rank. "
           "Defaults to 0 in most cases, when the main learner rank obtains all data. "
           "If run had MPI ranks collecting episodes it will default to 1.")
    parser.add_argument('--agentID', type=int, default=-1,
      help="Will plot the cumulative rewards obtained by this agent within each "
           "environment. Defaults to plotting all agents.")
    parser.add_argument('--ask',   dest='bAsk',   action='store_true',
    help="Set if script should ask before plotting each reward file (useful for multi-agent sims).")
    parser.set_defaults(bAsk=False)
    parsed = parser.parse_args()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    colorID = 0
    lines, fills = [], []

    for PATH in parsed.folders:

        if os.path.isdir(PATH): ensambleList = [PATH]
        else: ensambleList = glob.glob(PATH)

        if len(ensambleList) == 0: assert False, 'Directory %s not found'%PATH

        if parsed.agentID < 0: # plot all agents
            agentID = 0
            while 1: # loop over possible agents
                rank = getWorkerRankToPlot(parsed, agentID, ensambleList[0])
                if rank < 0:
                    # no data found for this agent on any rank,
                    # probably reached last agent
                    break
                fill, line = plotReward(parsed, ax, ensambleList, agentID, rank, colorID)
                agentID += 1
                if fill is not None:
                    colorID += 1
                    lines += [line]
                    fills += [fill]
        else:
            rank = getWorkerRankToPlot(parsed, parsed.agentID, ensambleList[0])
            fill, line = plotReward(parsed, ax, ensambleList, parsed.agentID, rank, colorID)
            if fill is not None:
                    colorID += 1
                    lines += [line]
                    fills += [fill]

    ax.set_xlabel('time steps')
    ax.set_ylabel('average cumulative rewards')
    ax.legend(loc='lower right')
    plt.tight_layout()
    plt.show()

  
