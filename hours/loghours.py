import os
from time import strftime

class Time(object):
    def __init__(self, years = None, months = None, days = None, hours = None, minutes = None, seconds = None):
        self.years = years if years is not None else strftime("%Y")
        self.months = months if months is not None else strftime("%m")
        self.days = days if days is not None else strftime("%d")
        self.hours = hours if hours is not None else strftime("%H")
        self.minutes = minutes if minutes is not None else strftime("%M")
        self.seconds = seconds if seconds is not None else strftime("%S")

        self._vals = [int(self.years), int(self.months), int(self.days), int(self.hours), int(self.minutes), int(self.seconds)]

    def get_formatted_time(self):
        return self.months + '/' + self.days + '/' + self.years + ' ' + self.hours + ':' + self.minutes + ':' + self.seconds

    def calculate_difference(self, other_time):
        diffs = []
        for i, j in zip(self._vals, other_time._vals):
            diffs.append(abs(i - j))
        


def main():
    ticks = Time()
    print(ticks.get_formatted_time())

if __name__ == '__main__':
    main()
