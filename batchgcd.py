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

def batchgcd(xs):
    def _product(xs, a, b):
        p = gmpy2.mpz(1)
        for k in range(a, min(len(xs), b)):
            p = gmpy2.mul(p, xs[k])
        return p

    def _product_tree(xs):
        tree = [xs]
        while len(xs) > 1:
            LOGGER.info('Calculating product tree: %10d' %(len(xs)))
            xs = [_product(xs, k * 2, (k + 1) * 2)
                for k in range((len(xs) + 1)//2)]
            tree.append(xs)
        return tree

    def _batchgcd(xs):
        tree = _product_tree(xs)
        rems = tree.pop()
        while tree:
            LOGGER.info('Calculating batch GCDs:   %10d' %(len(tree)))
            xs = tree.pop()
            rems = [gmpy2.mod(rems[i//2], gmpy2.mul(xs[i], xs[i]))
                    for i in range(len(xs))]
        return [gmpy2.gcd(gmpy2.t_div(r, n), n)
                for r, n in zip(rems, xs)]

    result = _batchgcd(xs)
    for i, x in enumerate(result):
        if x != xs[i]: # success
            yield x
        else: # fallback
            # only check where bgcd != 1 and x != y
            for d in (gmpy2.gcd(x, y) for j, y in enumerate(result)
                      if xs[j] != 1 and x != y):
                if d != 1:
                    yield d
                    break
            else:
                yield 1

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
