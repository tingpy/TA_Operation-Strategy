{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Problem 1"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Solves Wagner-Whitin problem using dynamic programming (DP).\n",
    "\n",
    "INPUT ARGUMENTS: \n",
    "   K       = fixed cost [scalar or vector]\n",
    "   h       = holding cost per item per period [scalar or vector]\n",
    "   T       = number of periods in horizon [scalar or vector]\n",
    "   d       = demand in each period [vector or matrix]\n",
    "   I0      = initial inventory [scalar or vector]\n",
    "\n",
    "OUTPUT ARGUMENTS:\n",
    "   Q       = optimal order quantities in each period [vector or matrix]\n",
    "   optcost = optimal cost [scalar or vector]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 65,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "WW_DP (generic function with 1 method)"
      ]
     },
     "execution_count": 65,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "function WW_DP(K, h, T, d, I0)\n",
    "    m = length(K)\n",
    "    \n",
    "    Q = zeros(m, maximum(T))\n",
    "    optcost = zeros(m, 1)\n",
    "    \n",
    "    for r in 1:m\n",
    "        theta = zeros(1, maximum(T)+1)\n",
    "        opts = zeros(Int64, 1, maximum(T))\n",
    "\n",
    "        tt = 1\n",
    "        while I0[r] > 0\n",
    "            amt = min(I0[r], d[r,tt])\n",
    "            \n",
    "            d[r,tt] = d[r,tt] - amt\n",
    "            I0[r] = I0[r] - amt\n",
    "            \n",
    "            tt = tt + 1\n",
    "        end\n",
    "        \n",
    "        for tt in T[r]:-1:1\n",
    "            theta[tt] = 1.0e300\n",
    "\n",
    "            for s in (tt+1):(T[r]+1)\n",
    "                tempcost = K[r]\n",
    "\n",
    "                for i in tt:(s-1)\n",
    "                    tempcost = tempcost + h[r] * (i - tt) * d[r,i]\n",
    "                end\n",
    "\n",
    "                tempcost = tempcost + theta[s]\n",
    "\n",
    "                if tempcost < theta[tt]\n",
    "                    theta[tt] = tempcost\n",
    "                    opts[tt] = s\n",
    "                end\n",
    "\n",
    "            end\n",
    "\n",
    "        end\n",
    "\n",
    "        tt = 1\n",
    "        while d[r,tt] == 0\n",
    "            tt = tt + 1\n",
    "        end\n",
    "       \n",
    "        optcost[r] = theta[tt]\n",
    "        while tt < (T[r]+1)\n",
    "            Q[r,tt] = sum( d[ r, tt:(opts[tt]-1) ] );\n",
    "            tt = opts[tt]\n",
    "        end\n",
    "\n",
    "    end\n",
    "    \n",
    "    return Q, optcost\n",
    "end"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 66,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "([210.0 0.0 … 0.0 0.0; 18.0 0.0 … 0.0 0.0; 0.0 0.0 … 850.0 1050.0], [1380.0; 183.0; 6600.0])"
      ]
     },
     "execution_count": 66,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "K = [500 150 1100]\n",
    "h = [2 1 15]\n",
    "T = [4 5 8]\n",
    "I0 = [0 3 2000]\n",
    "d = [90 120 80 70 0 0 0 0 \n",
    "     6 7 3 0 5 0 0 0 \n",
    "     550 1100 750 1400 400 850 850 1050]\n",
    "\n",
    "Q, optcost = WW_DP(K, h, T, d, I0)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 67,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "3×8 Matrix{Float64}:\n",
       " 210.0  0.0  150.0     0.0    0.0    0.0    0.0     0.0\n",
       "  18.0  0.0    0.0     0.0    0.0    0.0    0.0     0.0\n",
       "   0.0  0.0  400.0  1400.0  400.0  850.0  850.0  1050.0"
      ]
     },
     "execution_count": 67,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "Q"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 68,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "3×1 Matrix{Float64}:\n",
       " 1380.0\n",
       "  183.0\n",
       " 6600.0"
      ]
     },
     "execution_count": 68,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "optcost"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Problem 2"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "function WW_Perturb_Snyder(K, h, T, d, sigma, n)\n",
    "% Solves Wagner-Whitin problem for random perturbations of the data set\n",
    "% specified and plots the order quantities.\n",
    "%\n",
    "% INPUT ARGUMENTS: \n",
    "%   K       = fixed cost [scalar]\n",
    "%   h       = holding cost per item per period [scalar]\n",
    "%   T       = number of periods in horizon [scalar]\n",
    "%   d       = demand in each period [vector]\n",
    "%   sigma   = SD of random noise to be added to each demand value [scalar]\n",
    "%   n       = number of perturbations to test [scalar]\n",
    "%\n",
    "% OUTPUT ARGUMENTS:\n",
    "%   [none]\n",
    "    \n",
    "    % initialize internal arrays\n",
    "    Q = zeros(n,T);\n",
    "    d_perturbed = zeros(1,T);\n",
    "    \n",
    "    % loop through perturbations\n",
    "    for r = 1:n\n",
    "        \n",
    "        % perturb demands by adding N(0,sigma^2) random variate\n",
    "        d_perturbed = d + sigma * randn(1,T);\n",
    "        \n",
    "        % solve DP\n",
    "        [Q(r,:),temp] = WW_DP_Snyder(K, h, T, d_perturbed);\n",
    "        \n",
    "    end;\n",
    "    \n",
    "    % plot Q\n",
    "    figure;\n",
    "    plot(Q');\n",
    "        \n",
    "end"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "INPUT ARGUMENTS: \n",
    "\n",
    "   K       = fixed cost [scalar]\n",
    "   h       = holding cost per item per period [scalar]\n",
    "   T       = number of periods in horizon [scalar]\n",
    "   d       = demand in each period [vector]\n",
    "   sigma   = SD of random noise to be added to each demand value [scalar]\n",
    "   n       = number of perturbations to test [scalar]\n",
    "\n",
    "Output:\n",
    "\n",
    "   A figure"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "using PyPlot\n",
    "\n",
    "function WW_create_plot(K, h, T, d, sigma, n)\n",
    "    Q = zeros(n,T)\n",
    "    d_perturbed = zeros(1,T)\n",
    "    \n",
    "    for r in 1:n\n",
    "        \n",
    "        d_perturbed = d .+ sigma .* randn(1,T)\n",
    "        \n",
    "        % solve DP\n",
    "        [Q[r,:],temp] = WW_DP(K, h, T, d_perturbed)\n",
    "        \n",
    "    end\n",
    "    \n",
    "    figure\n",
    "    plot(Q')\n",
    "    \n",
    "end"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "fig = figure()\n",
    "\n",
    "# Plotting two datasets\n",
    "plot(iter, lower_bound, color=\"red\", linewidth=2.0, linestyle=\"-\",\n",
    " marker=\"o\", label=L\"Lower Bound $Z^k_L$\")\n",
    "plot(iter, upper_bound, color=\"blue\", linewidth=2.0, linestyle=\"-.\",\n",
    " marker=\"D\", label=L\"Upper Bound $Z^k_U$\")\n",
    "\n",
    "# Labeling axes\n",
    "xlabel(L\"iteration clock $k$\", fontsize=\"xx-large\")\n",
    "ylabel(\"objective function value\", fontsize=\"xx-large\")\n",
    "\n",
    "# Putting the legend and determining the location\n",
    "legend(loc=\"upper right\", fontsize=\"x-large\")\n",
    "\n",
    "# Add grid lines\n",
    "grid(color=\"#DDDDDD\", linestyle=\"-\", linewidth=1.0)\n",
    "tick_params(axis=\"both\", which=\"major\", labelsize=\"x-large\")\n",
    "\n",
    "# Title\n",
    "title(\"Lower and Upper Bounds\")\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Julia 1.6.1",
   "language": "julia",
   "name": "julia-1.6"
  },
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.6.1"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
