{{ fullname | escape | underline }}

.. automodule:: {{ fullname }}

{% block functions %}
{% if functions or classes %}
.. rubric:: Summary

.. autosummary::
   :nosignatures:
{% for item in classes %}
   {{ item }}
{%- endfor %}
{% for item in functions %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}

{% if classes %}
.. rubric:: Classes
{% for item in classes %}
.. autoclass:: {{ item }}
   :members:
{%- endfor %}
{% endif %}

{% if functions %}
.. rubric:: Functions
{% for item in functions %}
.. autofunction:: {{ item }}
{%- endfor %}
{% endif %}

{% block modules %}
{% if modules %}
.. rubric:: Modules

.. autosummary::
   :toctree:
   :recursive:
{% for item in modules %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}
