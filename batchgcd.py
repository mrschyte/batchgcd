import functools
import fileinput
import logging
import gmpy2

from multiprocessing import Pool, cpu_count

from cryptography import x509
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives.serialization import Encoding

from Crypto.PublicKey import RSA
from base64 import b64decode

logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger('batchgcd')

def calculate_r_mod_xsqr(params):
    r, x = params
    return r % x**2

def calculate_product(xs):
    return functools.reduce(lambda x, y: x * y, xs)

def batchgcd(xs):
    pool = Pool(processes=cpu_count())
    def _product_tree(xs):

        tree = [xs]
        while len(xs) > 1:
            LOGGER.info('Calculating product tree: %10d' %(len(xs)))
            #xs = [calculate_product(xs[k * 2 : (k + 1) * 2]) for k in range((len(xs) + 1)//2)]
            xs = pool.map(calculate_product, [xs[k * 2 : (k + 1) * 2] for k in range((len(xs) + 1)//2)]) 
            tree.append(xs)
        return tree

    def _batchgcd(xs):
        tree = _product_tree(xs)
        rems = tree.pop()
        while tree:
            LOGGER.info('Calculating batch GCDs:   %10d' %(len(tree)))
            xs = tree.pop()
            #rems = [rems[i//2] % xs[i]**2 for i in range(len(xs))]
            rems = pool.map(calculate_r_mod_xsqr, [(rems[i//2], xs[i]) for i in range(len(xs))])
        for r, n in zip(rems, xs):
            yield gmpy2.gcd(r//n, n)

    for k, x in enumerate(_batchgcd(xs)):
        if x != xs[k]: # success
            yield x
        else: # fallback
            for d in (gmpy2.gcd(x, y) for y in xs):
                if d != 1:
                    yield x // d
                    break

def load_single_certificate(data):
    try:
        cert = x509.load_der_x509_certificate(b64decode(data), backend=default_backend())
        pubk = cert.public_key().public_numbers()
        n, e = map(gmpy2.mpz, [pubk.n, pubk.e])
        return (data, n, e)
    except Exception as ex:
        LOGGER.debug(ex)
        return None

def load_certificates():
    moduli = []
    expns = []
    certs = []

    pool = Pool(cpu_count())
    temp = set()

    for result in pool.map(load_single_certificate, fileinput.input()):
        if result != None:
            data, n, e = result
            if n not in temp:
                certs.append(data)
                moduli.append(n)
                expns.append(e)
                temp.add(n)
    temp.clear()
    return (moduli, expns, certs)

def main():
    LOGGER.info('Loading certificates...')
    moduli, expns, certs = load_certificates()

    for k, q in enumerate(batchgcd(moduli)):
        if q != 1:
            n, e = int(moduli[k]), int(expns[k])
            p, q = int(n // q), int(q)
            d = int(gmpy2.invert(e, (p - 1) * (q - 1)))
            priv = RSA.construct((n, e, d))

            print("%s,%s" %(
                certs[k], ''.join(priv.exportKey().decode('ascii').split('\n')[1:-1])
            ))

if __name__ == "__main__":
    main()
