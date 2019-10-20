import matplotlib.pyplot as plt
import re
from pathlib import Path
class Iteration:
    def __init__(self, i):
        self.i = i 
        self.events = []
        return 
class event:
    def __init__(self, num, timestamp=None, local_rank=None, coordinate=None):
        self.num = num
        return 

def readfile(fpath='./output/1'):
    logfile = list(Path(fpath).glob('logfile*.txt'))[0]
    # timestamp_p = r'Timestamp : .*(\d{2}:\d{2}:\d{2}).*\n'

    iteration_p = 'Iteration (\d+$)\n'
    event_num_p = 'event number (\d+).*\n'

    log = []
    with logfile.open('r') as f:
        for i in f:
            # if re.match(timestamp_p, i):
            #     timestamp.append(re.search(timestamp_p, i).group(1))
            if re.match(iteration_p, i):
                t = int(re.search(iteration_p, i).group(1))
                new_iter = Iteration(t)
                log.append(new_iter)
            # print(event_num_p, i)
            if re.match(event_num_p, i):
                new_e = int(re.search(event_num_p, i).group(1))
                log[-1].events.append(event(new_e))
    return log
        #date = datetime.datetime.strptime(match.group(), '%Y-%m-%d').date()
log = readfile()
x = []
y = []
e_freq = []
for k in log:
    x.extend([k.i]*len(k.events))
    e_freq.append(len(k.events))
    for e in k.events:
        y.append(e.num)
# print(len(x), len(y))
# plt.hist(x)
# plt.plot(x, y)
# print(x)

plt.hist(y, density=False)
plt.xlabel('event')
plt.ylabel('total number of activation')

print(y)
plt.show()
# print(iteration)
# print(event_num)
# plt.plot(iteration, event_num)    

