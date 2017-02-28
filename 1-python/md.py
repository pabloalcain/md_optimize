import numpy as np
import pylab as pl
import itertools as it
import random as rn
import time

fp = open('time.dat', 'w')
t0 = time.clock()
timestep = 0.0005
n_particles = 100
n_steps = 1000
arr_potential = np.zeros(n_steps)
arr_kinetic = np.zeros(n_steps)
size = 6.0
position = np.zeros((n_particles, 3))
velocity = np.zeros((n_particles, 3))
force = np.zeros((n_particles, 3))
r_cutoff = 2.5
phi_cutoff = 4.0/r_cutoff**12 - 4.0/r_cutoff**6

number_side = int(np.ceil(np.cbrt(n_particles)))
distance = size/number_side
index_particle = 0

for i, j, k in it.product(range(number_side), range(number_side), range(number_side)):
  if index_particle == n_particles: break
  position[index_particle, 0] = i * distance
  position[index_particle, 1] = j * distance
  position[index_particle, 2] = k * distance
  index_particle += 1

for i in range(n_particles):
  for j in range(3):
    velocity[i, j] = rn.random()

potential = 0.0
force = np.zeros((n_particles, 3))
for i in range(n_particles):
  for j in range(n_particles):
    if i == j: continue
    dr = position[i, :] - position[j, :]
    for k in range(3):
      if dr[k] > size/2: dr[k] -= size
      elif dr[k] < -size/2: dr[k] += size
    distance = np.sqrt(np.dot(dr, dr))
    if distance < r_cutoff:
      force[i, :] += 24/distance**2 * (2.0/distance**12 - 1.0/distance**6) * dr
      potential += (4.0 * (1/distance**12 - 1/distance**6) - phi_cutoff)/2

kinetic = np.dot(velocity.flatten(), velocity.flatten())/2
print potential, kinetic, potential+kinetic

for l in xrange(n_steps):
  for i in range(n_particles):
    for j in range(3):
      if position[i, j] > size: position[i, j] -= size
      if position[i, j] < -size: position[i, j] += size
  position += velocity * timestep + force/2 * timestep**2
  velocity += force/2 * timestep

  potential = 0.0
  force = np.zeros((n_particles, 3))
  for i in range(n_particles):
    for j in range(n_particles):
      if i == j: continue
      dr = position[i, :] - position[j, :]
      for k in range(3):
        if dr[k] > size/2: dr[k] -= size
        elif dr[k] < -size/2: dr[k] += size
      distance = np.sqrt(np.dot(dr, dr))
      if distance < r_cutoff:
        force[i, :] += 24/distance**2 * (2.0/distance**12 - 1.0/distance**6) * dr
        potential += (4.0 * (1/distance**12 - 1/distance**6) - phi_cutoff)/2
  velocity += force/2 * timestep
  kinetic = np.dot(velocity.flatten(), velocity.flatten())/2
  arr_kinetic[l] = kinetic
  arr_potential[l] = potential
  print>>fp, l, time.clock() - t0
  print l, time.clock() - t0

pl.plot(arr_kinetic, label='Kinetic Energy')
pl.plot(arr_potential, label='Potential Energy')
pl.plot(arr_kinetic + arr_potential, label='Total Energy')
pl.savefig('energia.pdf')
pl.show()
