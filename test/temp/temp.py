# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

class A(object):     # deriving from 'object' declares A as a 'new-style-class'
    def foo(self):
        print("foo")

class B(A):
    def foo(self):
        super(B, self).foo()   # calls 'A.foo()'

myB = B()
myB.foo()