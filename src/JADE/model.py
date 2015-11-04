# __author__ = 'WeiFu'
"""
Test Functions here!
"""
import random
import math


class BaseModel(object):
  """
  base class
  """

  def __init__(self, hi=100, lo=-100):
    self.dimension = 30
    self.max = hi
    self.min = lo

  def generate(self):
    """
    generate new solution
    """
    return [random.uniform(self.min, self.max) for _ in xrange(self.dimension)]

  def evaluate(self, lst):
    """
    'abstract' evaluation method, the child class need to implement it
    """
    raise NotImplementedError("Evaluation!!")


class F1(BaseModel):
  """
  test F1
  """

  def __init__(self, hi=100, lo=-100):
    super(F1, self).__init__(hi, lo)

  def evaluate(self, lst):
    return sum([x ** 2 for x in lst])


class F2(BaseModel):
  """
  test F2
  """

  def __init__(self, hi=10, lo=-10):
    super(F2, self).__init__(hi, lo)

  def evaluate(self, lst):
    """
    evaluate the function
    """
    out = 1
    for one in lst:
      out *= math.fabs(one)
    return sum([math.fabs(x) for x in lst]) + out


class F4(BaseModel):
  """
  test F4
  """

  def __init__(self, hi=100, lo=-100):
    super(F4, self).__init__(hi, lo)

  def evaluate(self, lst):
    """
    evaluate the function
    """
    return max([math.fabs(each) for each in lst])


class F6(BaseModel):
  """
  test F6
  """

  def __init__(self, hi=100, lo=-100):
    super(F6, self).__init__(hi, lo)

  def evaluate(self, lst):
    """
    evaluate the function
    """
    return sum([math.floor(x + 0.5) ** 2 for x in lst])


class F7(object):
  """
  test F7
  """

  def __init__(self, hi=1.28, lo=-1.28):
    super(F7, self).__init__(hi, lo)

  def evaluate(self, lst):
    """
    evaluate the function
    """
    return sum([(i + 1) * (x) ** 4
                for i, x in enumerate(lst)]) + random.random()


class F9(object):
  """
  test F9
  """

  def __init__(self, hi=5.12, lo=-5.12):
    super(F9, self).__init__(hi, lo)

  def evaluate(self, lst):
    """
    evaluate the function
    """
    return sum([x * x - 10 * math.cos(2 * math.pi * x) + 10 for x in lst])
