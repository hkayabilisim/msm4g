==========
Data Types
==========

.. c:type:: LinkedListElement

   A node in the linked list deneme.
.. code-block:: c

    typedef struct LinkedListElement
    {
        void *data;
        struct LinkedListElement *next;
        struct LinkedListElement *prev;
    } LinkedListElement;

.. c:type:: LinkedList

   The main data structure used to store a collection
   of elements of type :c:type:`LinkedListElement`.
   It basically contains pointers to the head and tail elements:

.. code-block:: c

    typedef struct LinkedList
    {
     struct LinkedListElement *head;
     struct LinkedListElement *tail;
    } LinkedList;

.. autosummary::
   :nosignatures:

   sphinx.environment.BuildEnvironment
   sphinx.util.relative_uri