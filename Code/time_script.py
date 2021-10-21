import time
import sys

from multiprocessing import Process

integer=sys.argv[1]

init=map(int, integer.strip('[]'))

num=list(init)[0]

def exclaim(int):

    time.sleep(int)

    print("You were very patient!")

if __name__ == '__main__':

    program = Process(target=exclaim, args=(num,))

    program.start()

    program.join(timeout=5)

    program.terminate()
