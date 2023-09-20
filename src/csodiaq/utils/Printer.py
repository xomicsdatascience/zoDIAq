from timeit import default_timer as timer
import time


class Printer:
    _singletonInstance = None

    def __new__(self):
        if not self._singletonInstance:
            self._singletonInstance = super(Printer, self).__new__(self)
            self.startTime = timer()
            self.timeSinceLastPrint = 0
            self.timeSinceLastCheckpoint = 0
        return self._singletonInstance

    def __call__(self, printText, checkPoint=False):
        currentTime = timer()
        fiveMinutes = 300
        if checkPoint and self.timeSinceLastPrint < fiveMinutes:
            self.timeSinceLastPrint += (currentTime - self.timeSinceLastCheckpoint)
            self.timeSinceLastCheckpoint = currentTime
            return
        self.timeSinceLastCheckpoint = currentTime
        self.timeSinceLastPrint = 0
        readableCurrentTime = self.determine_time_since_start_in_human_readable_format(
            currentTime
        )
        print(
            "\n----> " + str(readableCurrentTime) + ": " + printText + "\n", flush=True
        )

    def determine_time_since_start_in_human_readable_format(self, currentTime):
        return time.strftime("%H:%M:%S", time.gmtime(currentTime - self.startTime))
