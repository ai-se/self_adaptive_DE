from __future__ import division
import random
import pdb
import numpy as np
# import math
from scipy.stats import cauchy
from model import *


class DeBase(object):
  def __init__(self, model, repeats=1500):
    self.model = model
    self.np = 100
    self.p = int(100 * 0.05)
    self.u_cr = 0.5  # initial mean of normal distribution
    self.u_f = 0.5  # initial mean of cauchy distribution
    self.gen_f()
    self.gen_cr()
    self.S_f = []  # list to store successful f's
    self.S_cr = []  # list to store successful cr's
    self.repeats = repeats
    self.obj = 1  # 1 is for minimization
    self.evaluation = 0
    self.scores = {}
    self.frontier = [model.generate() for _ in xrange(self.np)]
    self.evaluate()
    self.bestconf, self.bestscore = self.best()

  def evaluate(self):
    for n, arglst in enumerate(self.frontier):
      self.scores[n] = self.model.evaluate(arglst)

  def best(self):
    sortlst = []
    if self.obj == 1:  # this is for pf
      sortlst = sorted(self.scores.items(), key=lambda x: x[-1], reverse=True)  # alist of turple
    else:
      sortlst = sorted(self.scores.items(), key=lambda x: x[-1])  # alist of turple
    bestconf = [self.frontier[sortlst[i][0]] for i in xrange(-1 * self.p, 0)]
    bestscore = sortlst[-1][-1]
    return bestconf, bestscore

  def callModel(self, lst):
    raise NotImplementedError("callMode error")

  def treat(self, lst):
    """
    some parameters may have constraints, for example:
    when generating a parameter list, p[4]should be greater than p[5]
    You should implement this function in subclass
    """
    return NotImplementedError("treat error")

  def trim(self, x):
    return max(self.model.min, min(x, self.model.max))

  def gen3(self, n, f):
    seen = [n]

    def gen1(seen):
      while 1:
        k = random.randint(0, self.np - 1)
        if k not in seen:
          seen += [k]
          break
      return self.frontier[k]

    a = gen1(seen)
    b = gen1(seen)
    c = gen1(seen)
    return a, b, c

  def update(self, index, old):
    """
    generate new candidates
    """
    newf = []
    a, b, c = self.gen3(index, old)
    bestp = self.pickbest()
    for k in xrange(len(old)):
      newf.append(old[k] if self.cr[index] < random.random() else self.trim(
        (old[k] + self.fa[index] * (bestp[k] - old[k]) + self.fa[index] * (b[k] - c[k]))))
      # newf.append(old[k] if self.cr[index] < random.random() or random.randint(0, len(old) - 1) == k else self.trim(
      #   (old[k] + self.fa[index] * (bestp[k] - old[k]) + self.fa[index] * (b[k] - c[k]))))
    return self.treat(newf)

  def pickbest(self):
    """
    get the best one randomly from top p%
    """
    index = random.randint(0, len(self.bestconf) - 1)
    return self.bestconf[index]

  def update_mean1(self, old_val, mean_val):
    c = 0.1  ## this is set by JADE author
    val = (1 - c) * old_val + c * mean_val
    return val

  def lehmer_mean(self, SF):
    """
    this to calculate Lehmer mean of SF
    """
    if not len(SF):
      return 0
    return float(sum([f ** 2 for f in SF])) / sum(SF)

  def arithmetic_mean(self, Scr):
    if not len(Scr):
      return 0
    return float(sum(Scr)) / len(Scr)

  def update_mean(self):
    self.u_cr = self.update_mean1(self.u_cr, self.arithmetic_mean(self.S_cr))
    self.u_f = self.update_mean1(self.u_f, self.lehmer_mean(self.S_f))

  def gen_cr(self):
    """
    Cr in [0,1]
    """
    lst = []
    while len(lst) < self.np:
      lst.append(max(0, min(1, random.gauss(self.u_cr, 0.1))))
    # generate self.np cr's for each candidate from cauchy(0.5, 1) as initialization
    self.cr = lst[:]
    return

  def gen_f(self):
    """
    F in (0,1]
    """
    lst = []
    while len(lst) < self.np:
      temp = cauchy(self.u_f, 0.1).rvs()
      if temp >= 1:
        lst.append(1)
      elif temp <= 0:
        continue
      else:
        lst.append(temp)
    self.fa = lst[:]
    return

  def de(self):
    """
    main de body
    """

    def better(new, old):
      """
      better
      """
      return new < old if self.obj == 1 else new > old

    for _ in xrange(self.repeats):
      nextgeneration = []
      self.S_cr = []  ## clear before each generation
      self.S_f = []  ## clear before each generation
      for index, candidate in enumerate(self.frontier):
        new = self.update(index, candidate)
        newscore = self.callModel(new)
        self.evaluation += 1
        if better(newscore, self.scores[index]):
          nextgeneration.append(new)
          self.S_cr.append(self.cr[index])  # add cr and fa if new candidate is better
          self.S_f.append(self.fa[index])
          self.scores[index] = newscore
        else:
          nextgeneration.append(candidate)
      self.frontier = nextgeneration[:]
      newbestconf, newbestscore = self.best()
      if better(newbestscore, self.bestscore):
        self.bestscore = newbestscore
        self.bestconf = newbestconf[:]
      # update the means of gaussian distrobitopm amd caushy distribution
      self.update_mean()
      self.gen_cr()  # generate the new Cr for next generation
      self.gen_f()  # generate the new F for next generation
    print "final bestescore : " + str(self.bestscore)
    print "DONE !!!!"
    return self.bestscore


class JADE(DeBase):
  def __init__(self, model, repeats):
    super(JADE, self).__init__(model, repeats)

  def treat(self, lst):
    return lst

  def callModel(self, lst):
    return self.model.evaluate(lst)


if __name__ == "__main__":
  # test = [(1500,F1), (2000, F2 )(5000,F4), (100, F6), (1500,F6),(3000, F7),(1000, F9)]
  test = [(2000, F2)]
  for testcase in test:
    random.seed(1)
    result = []
    _repeats = testcase[0]
    _model = testcase[1]
    for _ in xrange(50):
      X = JADE(_model(), _repeats)
      result.append(X.de())
    print _model().__class__.__name__ + " mean is:" + str(np.mean(result))
    print _model().__class__.__name__ + " std is:" + str(np.std(result))
