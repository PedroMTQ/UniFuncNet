import os

from unifuncnet.Utils.util import find_path,UniFuncNet_FOLDER,SPLITTER

#External modules
from bs4 import BeautifulSoup
import time
from os.path import isdir


import requests
from urllib.parse import quote_plus
from random import choice,randrange
import re
import sys
from unifuncnet.Utils.util import xstr




class wait_for_display(object):
    def __init__(self, locator,wanted_atr,atr_value):
        self.locator = locator
        self.wanted_atr = wanted_atr
        self.atr_value = atr_value


    def __call__(self, driver):
        try:
            element = EC._find_element(driver, self.locator)
            print('waiting for display',element)
            return element.get_attribute(self.wanted_atr) == self.atr_value
        except StaleElementReferenceException:
            return False

class Web_Connector():
    def __init__(self,politeness_timer=10,
                 retry_timer=5,
                 try_limit=1,
                 test=False,
                 omit_error_messages=False,
                 timeout=20,
                 multiplier_max_timer=3):
        self.broken_link_file=f'{UniFuncNet_FOLDER}broken_links.txt'
        self.test=test
        self.timeout=timeout
        self.multiplier_max_timer=multiplier_max_timer
        self.politeness_timer = politeness_timer
        self.initial_politeness_timer = float(politeness_timer)
        self.retry_timer = retry_timer
        self.try_limit = try_limit
        self.request_time = time.time()
        self.omit_error_messages=omit_error_messages
        if self.test and not self.omit_error_messages: print('THIS FETCHER IS FOR TESTING, REMEMBER TO CHANGE IT!')


    def get_timeout(self):
        return self.timeout

    def is_test(self):
        return self.test

    def omit_errors(self):
        return self.omit_error_messages

    def get_request_time(self):
        return self.request_time


###################################
###############TIMER###############
###################################


    def get_politeness_timer(self):
        return self.politeness_timer

    def get_initial_politeness_timer(self):
        return self.initial_politeness_timer

    def get_retry_timer(self):
        return self.retry_timer

    def get_try_limit(self):
        return self.try_limit

    def set_politeness_timer(self,politeness_timer):
            self.politeness_timer=politeness_timer

    def get_randrange_retry_timer(self,denominator=4):
        allowed_range=[self.get_retry_timer()-self.get_retry_timer()/denominator,self.get_retry_timer()+self.get_retry_timer()/denominator]
        return round(randrange(int(allowed_range[0]*100),int(allowed_range[1]*100))/100,2)

    #This ensures that time between requests accompanies database latency
    def dynamic_politeness_timer(self,req_start,req_stop,c):
        """
        :param req_start: Time request was created
        :param req_stop: Time response was received
        :param c: Number of retries done so far
        :return: sets a politness timer based on the number of retries and on latency. Adding randomness for harder detection
        """
        if not req_start or not req_stop: return None
        try_limit=self.get_try_limit()
        ratio= int((c/try_limit)*100)
        multiplier=1
        if ratio in range(25,50):   multiplier=2
        elif ratio in range(50,75): multiplier=5
        if ratio > 75:              multiplier=10
        latency= req_stop-req_start
        if latency> self.get_politeness_timer():
            #pages may not give a status code caught by the request function. An unexistant page would then give a very high latency
            #we use latency<self.get_initial_politeness_timer()*3 as the threshold but it could be bigger
            if self.initial_politeness_timer < latency < self.initial_politeness_timer*self.multiplier_max_timer:
                rand_time=(randrange(0,100)/100) * latency + latency
                self.set_politeness_timer(rand_time*multiplier)
        else:
            rand_range =randrange(0,25)/100
            rand_choice=choice([rand_range,-rand_range])
            #in order to not go below the minimum time
            if latency>self.initial_politeness_timer:
                rand_time=rand_choice * latency + latency
                self.set_politeness_timer(rand_time*multiplier)
            else:
                rand_time =  rand_choice * self.initial_politeness_timer + self.initial_politeness_timer
                self.set_politeness_timer(rand_time*multiplier)
        #print('timer set for', self.get_politeness_timer())

    def set_retry_timer(self,retry_timer):
        self.retry_timer=retry_timer

    def set_try_limit(self,try_limit):
        self.try_limit=try_limit

    def update_request_time(self):
        self.request_time=time.time()

    def allow_request(self):
        #randrange for unpredictable request times
        if time.time()-self.get_politeness_timer()-self.get_politeness_timer()*(randrange(0,100)/100)>=self.get_request_time():
            self.update_request_time()
            return True
        else:
            return False

###################################
########EXCEPTION HANDLING#########
###################################

    #sometimes response 200 is returned , even though the page actually responds as a "bad request"
    def proper_response(self,req,url,c):
        if '400 Bad Request' in req:
            return False,url
        if 'The requested URL is malformed.' in req:
            return False,url
        if 'ncbi' in url and 'Error Blocked Diagnostic' in req:
            return False,url
        if 'kegg.jp' in url and 'No such data was found.' in req:
            url=url.replace('kegg.jp','genome.jp')
            return False,url
        if 'genome.jp' in url and 'No such data was found.' in req:
            url=url.replace('genome.jp','kegg.jp')
            return False,url
        if 'hmdb' in url and '<title>We\'re sorry, but something went wrong (500)</title>' in req:
            return False,url
        return True,url


    def is_broken_link(self,page,link,current_c):
        try:
            page_source=page.text
        except:
            page_source=page.page_source
        if link and page_source:
            broken_link_pattern = re.compile('(404 - File or directory not found.)|'
                                             '(Invalid URL parameters)|'
                                             '(Failed to retrieve sequence)|'
                                             '(File not found)|'
                                             '(PUGREST\.NotFound)')
            if re.search(broken_link_pattern, page_source):
                #just to avoid executing the save command over and over
                if current_c+5>=self.try_limit:
                    if not self.omit_errors(): print('Saved as invalid URL '+link)
                return True
        return False

    def is_api(self,url):
        if      'rest.kegg.jp' or\
                'eutils.ncbi.nlm.nih.gov' or\
                'www.uniprot.org'\
                in url: return True


###################################
#############REQUESTS##############
###################################

    def api_KEGG(self,to_search, database=None, api_type='find',kegg_option=None):
        # kegg compound database= cpd
        # kegg reaction database= rn
        # kegg gene database= gn
        if isinstance(to_search,list):
            to_search[0]=quote_plus(to_search[0])
            to_search=''.join(to_search)
        else: to_search = quote_plus(to_search)
        if database:
            if to_search:
                query = api_type + '/' + database + '/' + to_search
                if kegg_option:
                    query+='/'+kegg_option
            else:
                query = api_type + '/' + database
        else:
            query = api_type + '/' + to_search
        url = f'http://rest.kegg.jp/{query}'
        webpage = self.try_until_catch(url)
        if not webpage: return []
        soup = BeautifulSoup(webpage, 'lxml')
        if not soup.contents: return []  # when search doesn't return anything
        result = soup.body.p
        if result:
            result = result.text.split('\n')
        else:
            return []
        return result[:-1]


    def print_status_code(self,status_code):
        if status_code   == 400:      print('Bad request')
        elif status_code == 401:    print('Unauthorized')
        elif status_code == 403:    print('Forbidden')
        elif status_code == 404:    print('Not found')
        elif status_code == 408:    print('Request timeout')
        elif status_code == 429:    print('Too many requests')
        elif status_code == 500:    print('Internal Server Error')
        else:                       print('Client error response')


    def try_until_catch(self,url, data=None, original_response=False,exceptional_try_limit=None):
        req_start,req_end=None,None
        while not self.allow_request():
            time.sleep(0.1)
        url = url.replace(' ', quote_plus(' '))
        error = False
        c = 0
        if exceptional_try_limit: try_limit=exceptional_try_limit
        else: try_limit=self.try_limit
        while c <= try_limit:
            header = {'User-Agent':'Mozilla/5.0'}
            try:
                req_start=time.time()
                if data:
                    #data is sent as a dictionary
                    if self.get_try_limit()-c<= 5 or self.is_test() or self.is_api(url):
                        req = requests.post(url, headers=header, data=data,timeout=self.get_timeout())
                    else:
                        req = requests.post(url, headers=header, data=data,timeout=self.get_timeout())
                else:
                    if self.get_try_limit()-c<= 5 or self.is_test() or self.is_api(url):
                        req = requests.get(url, headers=header,timeout=self.get_timeout())
                    else:
                        req = requests.get(url, headers=header,timeout=self.get_timeout())
                proper_response,url=self.proper_response(req.text, url,c)
                if not proper_response:
                    c+=1
                    raise ConnectionError
                if self.is_broken_link(req,url,c): c+=5
                if error and req.status_code == 200:
                    if not self.omit_errors(): print(url+' was successfully retrieved after '+str(c)+' retries.')
                #to set a timer which follows current website latency
                req_end=time.time()
                self.dynamic_politeness_timer(req_start,req_end,c)
                #END GOAL
                if req.status_code== 200:
                    if original_response: return req
                    else:                 return req.text
                else:
                    c+=1
                    if req.status_code==429:
                        time.sleep(self.politeness_timer)
                    if not self.omit_errors():
                        self.print_status_code(req.status_code)
            #previous status are the most common, but there can be exceptions
            except:
                req_end=time.time()
                c+=1
                self.dynamic_politeness_timer(req_start,req_end,c)
                randrange_retry_timer=self.get_randrange_retry_timer()
                if not self.omit_errors(): print('Server error (requests-try '+str(c)+') while getting ' + url + ' , retrying again in ' + str(randrange_retry_timer) + ' seconds.')
                time.sleep(randrange_retry_timer)
                error = True
        print('Couldn\'t open link: '+url)
        with open(self.broken_link_file,'a+') as file:
            file.write(f'{url}\n')
        return None



if __name__ == '__main__':
    f=Web_Connector()
    url='https://biocyc.org/META/NEW-IMAGE?type=EC-NUMBER&object=EC-1.1.1.370'
    a=f.try_until_catch(url)
    print(a)
