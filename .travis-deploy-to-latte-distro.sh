#! /bin/sh
# To be run after "make dist" (or "make distcheck").  Upload the tar.gz to latte-distro by making a git commit.
set -e
echo "DISTNAME=@PACKAGE@-@VERSION@.tar.gz; DISTPATTERN=@PACKAGE@-*.tar.gz; VERSION=@VERSION@; prefix=@prefix@; srcdir=@srcdir@; abs_srcdir=@abs_srcdir@" > .dist-vars.in
./config.status --quiet --file=.dist-vars
. .dist-vars
rm -f .dist-vars.in .dist-vars
## UPLOAD to github repository latte-distro
# This uses a secret deploy key.
# The key pair was created using:
#   ssh-keygen -t dsa -f .travis_ci_latte_distro_deploy_key
# The public key was uploaded to https://github.com/latte-int/latte-distro/settings/keys
# The private key was encrypted using
#   travis encrypt-file --com .travis_ci_latte_distro_deploy_key
# and appears as .travis_ci_latte_distro_deploy_key.enc in the repository.
# This also uploaded decryption keys in the form of environment variables on Travis CI.
if test x"$encrypted_48054582d2ae_key" != x -a x"$encrypted_48054582d2ae_iv" != x -a ! -r $srcdir/.travis_ci_latte_distro_deploy_key ; then
    openssl aes-256-cbc -K $encrypted_48054582d2ae_key -iv $encrypted_48054582d2ae_iv -in .travis_ci_latte_distro_deploy_key.enc -out .travis_ci_latte_distro_deploy_key -d
fi
if test -r $srcdir/.travis_ci_latte_distro_deploy_key; then 
    echo "Deployment key exists, attempting to upload"
    chmod 0600 $srcdir/.travis_ci_latte_distro_deploy_key
    export GIT_SSH_COMMAND="ssh -i $abs_srcdir/.travis_ci_latte_distro_deploy_key -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no"
    rm -Rf latte-distro-deploy
    git clone --depth 1 git@github.com:latte-int/latte-distro.git latte-distro-deploy
    rm -f latte-distro-deploy/$DISTPATTERN
    cp $DISTNAME latte-distro-deploy/
    sed -i.orig 's/^LATTE_VERSION=[$]/#LATTE_VERSION=$/g;s/^LATTE_VERSION=[0-9]*.*/LATTE_VERSION='${VERSION}'/g' Makefile.am
    msg="Upgrade LattE to $VERSION."
    msg="${msg} (Automatic upload from Travis CI, ${TRAVIS_REPO_SLUG} job=${TRAVIS_JOB_NUMBER} branch=${TRAVIS_BRANCH}"
    if [[ -n ${TRAVIS_TAG} ]]; then
        msg="${msg} tag=${TRAVIS_TAG}"
    fi
    msg="${msg} commit=${TRAVIS_COMMIT}"
    if [[ -n ${TRAVIS_PULL_REQUEST} ]]; then
        msg="${msg} pull_request=${TRAVIS_PULL_REQUEST}"
    fi
    msg="${msg})"
    (cd latte-distro-deploy && git --version && git add --all && git commit -m "${msg}" && git push)
fi
set +e
